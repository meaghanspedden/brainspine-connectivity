% Plot subject-specific peaks with KDE 
subs = {'OP00212', 'OP00213',  'OP00215', 'OP00219', ...
        'OP00220', 'OP00221', 'OP00224', 'OP00225', 'OP00226'};
which_ori = 'all';
HFC=1;
bandw=15; %mm smoothing for kde

%% Step 1: Load data and identify peak locations
peak_y_coords = nan(length(subs), 1);
significant_peak_y_coords = nan(length(subs), 1);  % track significant peaks only

for i = 1:length(subs)
    sub = subs{i};
    save_dir = fullfile('D:\MSST001', [sub '_contrast']);
    if HFC
    savename = sprintf('spine_source_%s_%s.mat', sub, which_ori);
    else
    savename = sprintf('spine_source_%s_%s_nohfc.mat', sub, which_ori);

    end
    loaded = load(fullfile(save_dir, savename));

%         if any(loaded.pSig == 1)
%         fprintf('Subject %s has significant peaks parametric.\n', sub);
%         end

    y_coords = loaded.src.pos(:, 2);
    tvals = loaded.tvals;
    thresholds = loaded.thresholds;

    % Find peak t-value location (most negative t-value)
    [~, idx] = min(tvals);
    peak_y_coords(i) = y_coords(idx);

    % Check significance
    if ~isnan(tvals(idx)) && ~isnan(thresholds(idx)) && tvals(idx) < thresholds(idx)
        significant_peak_y_coords(i) = y_coords(idx);
    end
end

fprintf('%g subjects with significant peak\n', sum(~isnan(significant_peak_y_coords)))

%% Step 2: KDE and plot
sig_mask = ~isnan(significant_peak_y_coords);

% KDE for all peaks
[f_all, xi_all] = ksdensity(peak_y_coords, 'Bandwidth', bandw);

% Generate unique colors for each subject
color_map = lines(length(subs));  % built-in colormap with distinct colors

% Create figure
figure;
plot_handle_all = plot(xi_all, f_all, 'color', [0.6, 0.6, 0.6], 'LineWidth', 2, 'LineStyle','--');
hold on;

% Plot individual peaks and store handles for legend
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

% KDE for significant peaks (if any)
if any(sig_mask)
    [f_sig, xi_sig] = ksdensity(significant_peak_y_coords(sig_mask), 'Bandwidth', bandw);
    plot_handle_sig = plot(xi_sig, f_sig, 'color', [0.24, 0.70, 0.44], 'LineWidth', 2);  % green KDE for sig

    % Combine all legend entries
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
title(sprintf('Peak negative t-values ori %s %s mm smoothing',which_ori,bandw));



successes=sum(~isnan(significant_peak_y_coords));

p = binocdf(successes-1, length(subs), 0.05, 'upper');  % one-tailed

fprintf('binomial p=%g\n',p)

%% proportion significant results greater than chance
