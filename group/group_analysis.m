subs = {'OP00212', 'OP00213', 'OP00214', 'OP00215', 'OP00219', ...
        'OP00220', 'OP00221', 'OP00224', 'OP00225', 'OP00226'};
which_ori = 'all';

%% Step 1: Load data and identify peak locations (no interpolation)
peak_y_coords = nan(length(subs), 1);
significant_peak_y_coords = nan(length(subs), 1);  % track significant peaks only

for i = 1:length(subs)
    sub = subs{i};
    save_dir = fullfile('D:\MSST001', [sub '_contrast']);
    savename = sprintf('spine_source_%s_%s.mat', sub, which_ori);

    loaded = load(fullfile(save_dir, savename));

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

%% Step 2: KDE and plot
sig_mask = ~isnan(significant_peak_y_coords);

% KDE for all peaks
[f_all, xi_all] = ksdensity(peak_y_coords, 'Bandwidth', 3);

% Generate unique colors for each subject
color_map = lines(length(subs));  % built-in colormap with distinct colors

% Plot all peaks KDE
figure;
plot(xi_all, f_all, 'color', [0.6, 0.6, 0.6], 'LineWidth', 2, 'LineStyle','--'); hold on;

% Plot individual peaks: color-coded by subject, transparency = significance
for i = 1:length(subs)
    y = peak_y_coords(i);
    is_sig = ~isnan(significant_peak_y_coords(i));
    color = color_map(i, :);
    alpha_val = 1.0 * is_sig + 0.4 * ~is_sig;  % 1 if significant, 0.4 otherwise

    scatter(y, 0, 80, 'v', 'filled', ...
        'MarkerFaceColor', color, ...
        'MarkerEdgeColor', color, ...
        'MarkerFaceAlpha', alpha_val, ...
        'MarkerEdgeAlpha', alpha_val);
end

% KDE for significant peaks (if any)
if any(sig_mask)
    [f_sig, xi_sig] = ksdensity(significant_peak_y_coords(sig_mask), 'Bandwidth', 3);
    plot(xi_sig, f_sig, 'color', [0.24, 0.70, 0.44], 'LineWidth', 2);  % green KDE for sig
    legend({'All Peaks (KDE)', 'Significant Peaks (KDE)'}, 'Location', 'best');
else
    warning('No significant peaks found. Only plotting all peaks.');
    legend({'All Peaks (KDE)'}, 'Location', 'best');
end

xlabel('Spinal Cord Position (mm)');
ylabel('Density');
title('Peak negative t-values (color-coded by subject)');
box off;
