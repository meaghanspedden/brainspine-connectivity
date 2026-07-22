%% compare_bem_bs_sensitivity_all_ori.m
% Computes R² sensitivity for both BEM and Biot-Savart leadfields
% across sensor-shifted conditions, for all 3 dipole orientations.

geoms_path  = 'C:\Users\mspedden\Documents\for_meaghan\geom_sens_shift';
bem_lf_path = 'C:\Users\mspedden\Documents\for_meaghan\bem_fields';
bs_lf_path  = 'C:\Users\mspedden\Documents\bslaw_sensitivity_analysis\bslaw_sensitivity_analysis\bs_law_fields';

%% ── Load BEM original ────────────────────────────────────────────────────────
orig_bem_dir  = fullfile(bem_lf_path, 'geometries_sensor_original');
orig_bem_file = dir(fullfile(orig_bem_dir, 'leadfield_*.mat'));
assert(~isempty(orig_bem_file), 'Could not find BEM original leadfield');

orig_bem    = load(fullfile(orig_bem_dir, orig_bem_file(1).name));
LF_bem_orig = extractLF_bem(orig_bem.leadfield_cord);
src_pos_bem = orig_bem.leadfield_cord.pos(orig_bem.leadfield_cord.inside, :);
cord_pos    = src_pos_bem(:, 2) * 1000;   % m → mm, inferior-superior axis
n_sources   = size(LF_bem_orig, 2);

%% ── Load BS original ─────────────────────────────────────────────────────────
orig_bs_file = fullfile(bs_lf_path, 'leadfield_sensor_original_bslaw_experimental.mat');
assert(exist(orig_bs_file, 'file') == 2, 'Could not find BS original leadfield');

orig_bs    = load(orig_bs_file);
LF_bs_orig = extractLF_bs(orig_bs.leadfield_bs);
src_pos_bs = orig_bs.leadfield_bs.pos;   % [53×3] mm

assert(size(LF_bs_orig, 2) == n_sources, ...
    'Source count mismatch: BEM=%d, BS=%d', n_sources, size(LF_bs_orig,2));

%% ── Discover shifted BEM folders ─────────────────────────────────────────────
bem_dirs = dir(fullfile(bem_lf_path, 'geometries_sensor_*'));
bem_dirs = bem_dirs([bem_dirs.isdir]);
bem_dirs = bem_dirs(~strcmp({bem_dirs.name}, 'geometries_sensor_original'));
fprintf('Found %d BEM sensor-shifted conditions\n', numel(bem_dirs));

%% ── Loop over shifted conditions ─────────────────────────────────────────────
results = struct();
row = 1;

for ci = 1:numel(bem_dirs)
    folder = bem_dirs(ci).name;
    label  = regexprep(folder, '^geometries_', '');

    bem_files = dir(fullfile(bem_lf_path, folder, 'leadfield_*_bem_experimental.mat'));
    if isempty(bem_files)
        warning('No BEM leadfield in %s — skipping', folder); continue
    end

    geom_file = fullfile(geoms_path, [folder '.mat']);
    if ~exist(geom_file, 'file')
        warning('No geometry file for %s — skipping', folder); continue
    end
    geom      = load(geom_file, 'shift_vec_mm');
    shift_vec = geom.shift_vec_mm;
    shift_mag = round(norm(shift_vec), 2);

    bs_file = fullfile(bs_lf_path, sprintf('leadfield_%s_bslaw_experimental.mat', label));
    if ~exist(bs_file, 'file')
        warning('No BS leadfield for %s — skipping', label); continue
    end

    fprintf('Processing: %s  |shift|=%.2f mm  (%d/%d)\n', ...
        label, shift_mag, ci, numel(bem_dirs));

    shifted_bem  = load(fullfile(bem_lf_path, folder, bem_files(1).name));
    LF_bem_shift = extractLF_bem(shifted_bem.leadfield_cord);

    shifted_bs  = load(bs_file);
    LF_bs_shift = extractLF_bs(shifted_bs.leadfield_bs);

    % R² for all 3 orientations
    r2_bem = zeros(n_sources, 3);
    r2_bs  = zeros(n_sources, 3);
    for ori = 1:3
        r2_bem(:, ori) = computeR2(LF_bem_orig(:,:,ori), LF_bem_shift(:,:,ori), n_sources);
        r2_bs(:, ori)  = computeR2(LF_bs_orig(:,:,ori),  LF_bs_shift(:,:,ori),  n_sources);
    end

    results(row).label      = label;
    results(row).shift_vec  = shift_vec;
    results(row).shift_mag  = shift_mag;
    results(row).r2_bem     = r2_bem;
    results(row).r2_bs      = r2_bs;

    fprintf('  BEM median R² — RL:%.4f  AP:%.4f  SI:%.4f\n', ...
        median(r2_bem(:,1)), median(r2_bem(:,2)), median(r2_bem(:,3)));
    fprintf('   BS median R² — RL:%.4f  AP:%.4f  SI:%.4f\n', ...
        median(r2_bs(:,1)),  median(r2_bs(:,2)),  median(r2_bs(:,3)));
    row = row + 1;
end

save(fullfile(bs_lf_path, 'sensitivity_comparison_BEM_vs_BS.mat'), ...
    'results', 'src_pos_bem', 'cord_pos');
fprintf('\nDone. %d conditions processed.\n', row-1);

%% ── Sensitivity summary at observed variability range (~15 mm) ───────────────
all_mags  = [results.shift_mag];
idx_15    = all_mags >= 12 & all_mags <= 18;
r2_bem_15 = vertcat(results(idx_15).r2_bem);
r2_bs_15  = vertcat(results(idx_15).r2_bs);

fprintf('\n── Sensitivity at observed variability range (12–18 mm) ──\n');
ori_labels = {'RL', 'AP', 'SI'};
for ori = 1:3
    fprintf('%s — BEM: median=%.4f  min=%.4f | BS: median=%.4f  min=%.4f\n', ...
        ori_labels{ori}, ...
        median(r2_bem_15(:,ori)), min(r2_bem_15(:,ori)), ...
        median(r2_bs_15(:,ori)),  min(r2_bs_15(:,ori)));
end

%% ── Figures: one per orientation ─────────────────────────────────────────────
uniq_mags = unique(all_mags);
n_mags    = numel(uniq_mags);

cmap_mags = interp1([0 1], [0.7 0.85 1.0; 0.02 0.15 0.50], ...
    linspace(0, 1, n_mags));

for ori = 1:3
    figure('Color','w','Position',[100 100 1200 480]);

    for panel = 1:2
        ax = subplot(1,2,panel); hold on;

        for mi = 1:n_mags
            idx_mag = find(all_mags == uniq_mags(mi));
            for ki = 1:numel(idx_mag)
                r = results(idx_mag(ki));
                if panel == 1
                    r2_vec = r.r2_bs(:, ori);   % BS on left
                else
                    r2_vec = r.r2_bem(:, ori);  % BEM on right
                end
                plot(cord_pos, r2_vec, ...
                    'Color', [cmap_mags(mi,:) 0.7], ...
                    'LineWidth', 1.5, ...
                    'HandleVisibility', 'off');
            end
        end

        yline(0.99, 'k:', 'LineWidth', 1.0, 'HandleVisibility', 'off');
        yline(0.95, 'k--', 'LineWidth', 1.0, 'HandleVisibility', 'off');
        text(min(cord_pos), 0.99, 'R²=0.99 ', 'FontSize', 14, ...
            'Color', [0.3 0.3 0.3], 'VerticalAlignment', 'bottom', ...
            'HorizontalAlignment', 'right');
        text(min(cord_pos), 0.95, 'R²=0.95 ', 'FontSize', 14, ...
            'Color', [0.3 0.3 0.3], 'VerticalAlignment', 'bottom', ...
            'HorizontalAlignment', 'right');

        ylim([0.8 1.02]);
        pad = 0.05 * range(cord_pos);
        xlim([min(cord_pos)-pad  max(cord_pos)+pad]);
        xlabel('Position along cord (mm, inferior → superior)', 'FontSize', 14);
        ylabel('R²', 'FontSize', 14);
        ax.FontSize = 14;
        panel_titles = {'Biot-Savart', 'BEM'};
        title(panel_titles{panel}, 'FontWeight', 'normal', 'FontSize', 14);
        grid on; box on;

        colormap(ax, interp1([0 1], [0.7 0.85 1.0; 0.02 0.15 0.50], linspace(0,1,256)));
        cb            = colorbar(ax);
        clim(ax, [min(uniq_mags) max(uniq_mags)]);
        cb.Label.String   = 'Shift magnitude (mm)';
        cb.Label.FontSize = 14;
        cb.FontSize       = 14;
    end

    sgtitle(sprintf('Sensor shift sensitivity — %s dipole orientation', ori_labels{ori}), ...
        'FontSize', 14);

    fig_out = fullfile(bs_lf_path, ...
        sprintf('sensitivity_comparison_BEM_vs_BS_%s.png', ori_labels{ori}));
    print(gcf, fig_out, '-dpng', '-r300');
    fprintf('Figure saved to %s\n', fig_out);
end

%% ── Simple 3-line summary figure ─────────────────────────────────────────────
%  Picks three representative shift magnitudes (small / medium / large),
%  averages R² across all conditions at that magnitude, and plots one line
%  per shift level for BS (left) and BEM (right). One figure per orientation.
%
%  Adjust target_shifts to match magnitudes present in your data.
target_shifts = [5, 15, 25];   % mm — nominal targets used only to SELECT a
                               % representative condition near each level
shift_cols    = [0.20 0.52 0.78;   % blue   — small
                 0.89 0.47 0.10;   % orange — medium
                 0.17 0.63 0.17];  % green  — large
shift_styles  = {'-', '--', ':'};  % solid / dashed / dotted

% Each target maps to the nearest AVAILABLE shift condition. Label the legend
% with that condition's TRUE magnitude, not the nominal target: the dataset
% contains no exact 5/15/25 mm shifts (largest available is ~20 mm), so a
% "25 mm" label would misstate the data.
actual_mags = zeros(size(target_shifts));
for ti = 1:numel(target_shifts)
    [~, ni]         = min(abs(all_mags - target_shifts(ti)));
    actual_mags(ti) = all_mags(ni);
end

for ori = 1:3
    figure('Color','w','Position',[100 100 1000 420]);

    if ori == 1
        fprintf('\n-- Simple figure: shift bundle summary --\n');
        fprintf('  %-15s  %-15s  %-10s  %s\n', ...
            'Target (mm)', 'Nearest (mm)', 'N conditions', 'Labels');
        for ti = 1:numel(target_shifts)
            tgt = target_shifts(ti);
            [~, nearest_idx] = min(abs(all_mags - tgt));
            nearest_mag      = all_mags(nearest_idx);
            idx_tgt          = find(all_mags == nearest_mag);
            lbls = strjoin({results(idx_tgt).label}, ', ');
            fprintf('  %-15d  %-15.2f  %-10d  %s\n', ...
                tgt, nearest_mag, numel(idx_tgt), lbls);
        end
        fprintf('\n');
    end

    for panel = 1:2
        ax = subplot(1,2,panel); hold on;

        leg_handles = gobjects(numel(target_shifts), 1);

        for ti = 1:numel(target_shifts)
            tgt = target_shifts(ti);

            % Find conditions closest to this target magnitude
            [~, nearest_idx] = min(abs(all_mags - tgt));
            nearest_mag      = all_mags(nearest_idx);
            idx_tgt          = find(all_mags == nearest_mag);

            % Min/max range across all conditions at this magnitude
            r2_mat = zeros(n_sources, numel(idx_tgt));
            for ki = 1:numel(idx_tgt)
                if panel == 1
                    r2_mat(:, ki) = results(idx_tgt(ki)).r2_bs(:, ori);
                else
                    r2_mat(:, ki) = results(idx_tgt(ki)).r2_bem(:, ori);
                end
            end
            r2_mean = mean(r2_mat, 2);
            r2_min  = min(r2_mat, [], 2);
            r2_max  = max(r2_mat, [], 2);

            % Shaded range
            patch([cord_pos; flipud(cord_pos)], ...
                  [r2_min; flipud(r2_max)], ...
                  shift_cols(ti,:), ...
                  'FaceAlpha', 0.15, 'EdgeColor', 'none', ...
                  'HandleVisibility', 'off');

            % Mean line on top
            leg_handles(ti) = plot(cord_pos, r2_mean, shift_styles{ti}, ...
                'Color',     shift_cols(ti,:), ...
                'LineWidth', 2.2);
        end

        yline(0.99, ':', 'Color',[0.3 0.3 0.3], 'LineWidth', 1.0, 'HandleVisibility','off');
        yline(0.95, '--','Color',[0.3 0.3 0.3], 'LineWidth', 1.0, 'HandleVisibility','off');

        ylim([0.8 1.02]);
        pad = 0.05 * range(cord_pos);
        xlim([min(cord_pos)-pad  max(cord_pos)+pad]);
        xlabel('Position along cord (mm, inferior → superior)', 'FontSize', 13);
        ylabel('R²', 'FontSize', 13);
        ax.FontSize = 13;
        panel_titles = {'Biot-Savart', 'BEM'};
        title(panel_titles{panel}, 'FontWeight','normal', 'FontSize', 13);
        grid on; box on;

        legend(leg_handles, ...
            arrayfun(@(m) sprintf('%.1f mm', m), actual_mags, 'UniformOutput', false), ...
            'Location','southwest', 'Box','off', 'FontSize', 11);
    end

    sgtitle(sprintf('Sensor shift sensitivity — %s dipole orientation', ori_labels{ori}), ...
        'FontSize', 13);

    fig_base = fullfile(bs_lf_path, ...
        sprintf('sensitivity_simple_BEM_vs_BS_%s', ori_labels{ori}));
    print(gcf, [fig_base '.png'], '-dpng', '-r300');
    exportgraphics(gcf, [fig_base '.pdf'], 'ContentType','vector');
    print(gcf, [fig_base '.svg'], '-dsvg', '-painters');
    fprintf('Simple figure saved: %s (.png/.pdf/.svg)\n', fig_base);
end

%% ── BS only, AP orientation, simplified 3-line figure ────────────────────────
ori_ap = 2;   % AP orientation index

figure('Color','w','Position',[100 100 520 420]);
hold on;

leg_handles_ap = gobjects(numel(target_shifts), 1);

for ti = 1:numel(target_shifts)
    tgt = target_shifts(ti);

    [~, nearest_idx] = min(abs(all_mags - tgt));
    nearest_mag      = all_mags(nearest_idx);
    idx_tgt          = find(all_mags == nearest_mag);

    r2_mat = zeros(n_sources, numel(idx_tgt));
    for ki = 1:numel(idx_tgt)
        r2_mat(:, ki) = results(idx_tgt(ki)).r2_bs(:, ori_ap);
    end
    r2_mean = mean(r2_mat, 2);
    r2_min  = min(r2_mat, [], 2);
    r2_max  = max(r2_mat, [], 2);

    patch([cord_pos; flipud(cord_pos)], ...
          [r2_min; flipud(r2_max)], ...
          shift_cols(ti,:), ...
          'FaceAlpha', 0.15, 'EdgeColor', 'none', ...
          'HandleVisibility', 'off');

    leg_handles_ap(ti) = plot(cord_pos, r2_mean, shift_styles{ti}, ...
        'Color', shift_cols(ti,:), 'LineWidth', 2.2);
end

yline(0.99, ':', 'Color',[0.3 0.3 0.3], 'LineWidth', 1.0, 'HandleVisibility','off');
yline(0.95, '--','Color',[0.3 0.3 0.3], 'LineWidth', 1.0, 'HandleVisibility','off');

ylim([0.8 1.02]);
pad = 0.05 * range(cord_pos);
xlim([min(cord_pos)-pad  max(cord_pos)+pad]);
xlabel('Position along cord (mm, inferior \rightarrow superior)', 'FontSize', 13);
ylabel('R²', 'FontSize', 13);
title('Biot-Savart — AP orientation', 'FontWeight','normal', 'FontSize', 13);
grid on; box on;
set(gca, 'FontSize', 13);

legend(leg_handles_ap, ...
    arrayfun(@(m) sprintf('%.1f mm', m), actual_mags, 'UniformOutput', false), ...
    'Location','southwest', 'Box','off', 'FontSize', 11);

fig_out = fullfile(bs_lf_path, 'sensitivity_bs_AP_simple.png');
print(gcf, fig_out, '-dpng', '-r300');
fprintf('BS AP simple figure saved: %s\n', fig_out);
function LF = extractLF_bem(leadfield_cord)
    idx       = find(leadfield_cord.inside);
    n_sources = numel(idx);
    n_sensors = numel(leadfield_cord.label);
    LF        = zeros(n_sensors, n_sources, 3);
    for k = 1:n_sources
        LF(:, k, :) = leadfield_cord.leadfield{idx(k)};
    end
end

function LF = extractLF_bs(leadfield_bs)
    n_sources = size(leadfield_bs.pos, 1);
    n_sensors = numel(leadfield_bs.label);
    LF        = zeros(n_sensors, n_sources, 3);
    for k = 1:n_sources
        LF(:, k, :) = leadfield_bs.leadfield{k};
    end
end

function r2 = computeR2(LF_orig, LF_shift, n_sources)
    r2 = zeros(n_sources, 1);
    for k = 1:n_sources
        r     = corrcoef(LF_orig(:,k), LF_shift(:,k));
        r2(k) = r(1,2)^2;
    end
end