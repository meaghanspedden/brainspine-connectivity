% Peak analysis across orientations: consistency and KDE

clear all
close all
subs = {'OP00212', 'OP00213',  'OP00215', 'OP00219', ...
    'OP00220', 'OP00221', 'OP00224', 'OP00225', 'OP00226'};

ori2test = {'1', '2', '3', 'all'};
HFC = 1;
bandw = 15; % mm smoothing for KDE

% Step 1: Initialize matrices to store peak y-coordinates
peak_y_coords_all = nan(length(subs), length(ori2test));
significant_peak_y_coords_all = nan(length(subs), length(ori2test));

for i = 1:length(subs)
    sub = subs{i};
    for oo = 1:length(ori2test)
        which_ori = ori2test{oo};
        save_dir = fullfile('D:\MSST001', [sub '_contrast']);
        
        if HFC
            savename = sprintf('spine_source_%s_%s.mat', sub, which_ori);
        else
            savename = sprintf('spine_source_%s_%s_nohfc.mat', sub, which_ori);
        end

        loaded = load(fullfile(save_dir, savename));
        y_coords = loaded.src.pos(:, 2);
        tvals = loaded.tvals;
        thresholds = loaded.thresholds;

        [~, idx] = min(tvals);
        peak_y_coords_all(i, oo) = y_coords(idx);

        if ~isnan(tvals(idx)) && ~isnan(thresholds(idx)) && tvals(idx) < thresholds(idx)
            significant_peak_y_coords_all(i, oo) = y_coords(idx);
        end
    end
end

% Step 2: Consistency metrics across all 4 orientations
peak_std = std(peak_y_coords_all, 0, 2, 'omitnan');    % std across all orientations
peak_range = range(peak_y_coords_all, 2);              % range across all orientations

fprintf('\nSubject-wise consistency (SD and range of peak y-coords across ori 1–2–3–all):\n');
for i = 1:length(subs)
    fprintf('%s: SD = %.2f mm, Range = %.2f mm\n', subs{i}, peak_std(i), peak_range(i));
end

% Optional: Threshold-based count
consistent_thresh = 10;
n_consistent = sum(peak_std < consistent_thresh);
fprintf('%d out of %d subjects had consistent peak locations (SD < %.1f mm)\n', ...
    n_consistent, length(subs), consistent_thresh);

% Step 3: KDE and plotting for 'all' orientation only
which_ori = 'all';
peak_y_coords = peak_y_coords_all(:, strcmp(ori2test, 'all'));
significant_peak_y_coords = significant_peak_y_coords_all(:, strcmp(ori2test, 'all'));

sig_mask = ~isnan(significant_peak_y_coords);
[f_all, xi_all] = ksdensity(peak_y_coords, 'Bandwidth', bandw);

% Generate unique colors for subjects
color_map = lines(length(subs));

% Create KDE plot
figure;
plot_handle_all = plot(xi_all, f_all, 'color', [0.6, 0.6, 0.6], 'LineWidth', 2, 'LineStyle','--');
hold on;

% Plot individual subject peaks
scatter_handles = gobjects(length(subs), 1);
for i = 1:length(subs)
    y = peak_y_coords(i);
    is_sig = ~isnan(significant_peak_y_coords(i));
    color = color_map(i, :);
    alpha_val = 1.0 * is_sig + 0.4 * ~is_sig;

    h = scatter(y, 0, 80, 'v', 'filled', ...
        'MarkerFaceColor', color, ...
        'MarkerEdgeColor', color, ...
        'MarkerFaceAlpha', alpha_val, ...
        'MarkerEdgeAlpha', alpha_val);
    
    scatter_handles(i) = h;
end

% KDE for significant peaks
if any(sig_mask)
    [f_sig, xi_sig] = ksdensity(significant_peak_y_coords(sig_mask), 'Bandwidth', bandw);
    plot_handle_sig = plot(xi_sig, f_sig, 'color', [0.24, 0.70, 0.44], 'LineWidth', 2);  % green

    legend([plot_handle_all, plot_handle_sig, scatter_handles'], ...
        ['All Peaks (KDE)', 'Significant Peaks (KDE)', subs], ...
        'Location', 'eastoutside');
else
    warning('No significant peaks found. Only plotting all peaks.');
    legend([plot_handle_all, scatter_handles'], ...
        ['All Peaks (KDE)', subs], ...
        'Location', 'eastoutside');
end

xlabel('Spinal Cord Position (mm)');
ylabel('Density');
title(sprintf('Peak negative t-values (%s orientation, %d mm smoothing)', which_ori, bandw));

% Step 4: Binomial test for significant results
successes = sum(~isnan(significant_peak_y_coords));
p = binocdf(successes-1, length(subs), 0.05, 'upper');  % one-tailed
fprintf('Binomial p = %g for %d/%d significant peaks\n', p, successes, length(subs));

% Step 5: Summary plots for consistency
figure;
boxplot(peak_std, 'Labels', {'Across 1,2,3,all'});
ylabel('Standard Deviation of Peak Y-Coords (mm)');
title('Consistency of Peak Locations Across Orientations');

% Step 6: Boxplot of peak variability per subject across orientations with stars for significance

figure;
boxplot(peak_y_coords_all', 'Labels', subs);  % transpose to get one box per subject
xlabel('Subject');
ylabel('Peak Y-Coordinate (mm)');
title('Peak Location Variability Across Orientations (per Subject)');
grid on;
hold on;

% Parameters for star plotting
numSubs = length(subs);
numOri = length(ori2test);

% Get y-axis limits to place stars a bit above boxes
ylims = ylim;
ystar = ylims(2) + 0.05 * range(ylims);  % 5% above top of y-axis

% We plot stars above each box: 
% boxplot groups are by subject, so x-axis is 1:numSubs
for subjIdx = 1:numSubs
    for oriIdx = 1:numOri
        if ~isnan(significant_peak_y_coords_all(subjIdx, oriIdx))
            % Plot star above the box for this subject
            % x-position is subject number
            % Add small horizontal offset for orientation so stars don't overlap
            
            % Offset: spread stars between -0.15 to +0.15 around the box center
            offset = ((oriIdx - (numOri+1)/2) / (numOri)) * 0.3;
            xPos = subjIdx + offset;
            
            % Plot star
            text(xPos, ystar, '*', 'FontSize', 16, 'HorizontalAlignment', 'center', 'Color', 'r');
        end
    end
end

% Adjust y limits to make room for stars
ylim([ylims(1), ystar + 0.05*range(ylims)]);