%% plot_spine_lf_topoplots.m
clear all; close all; clc;

%% =========================================================================
%  USER CONFIG
%% =========================================================================
fieldtrip_path = 'C:\Users\mspedden\Documents\fieldtrip';
spm_path       = 'C:\Users\mspedden\Documents\spm';
save_dir       = 'C:\Users\mspedden\Documents\brainspine_savetest\topoplots';

lf_bem_path = 'C:\Leadfields meshes\leadfield_experimental_bem_experimental.mat';
lf_bs_path  = 'C:\Leadfields meshes\leadfield_experimental_bslaw_experimental.mat';
geomfile    = 'C:\Leadfields meshes\geometries_experimental.mat';

source_idx = 15;   % for topoplots

% Orientations to analyse — col 1 = Left-Right (X), col 3 = Ventral-Dorsal (Z/radial)
ori_cols  = [1, 3];
ori_names = {'Left-Right', 'Ventral-Dorsal'};

%% =========================================================================
%  SETUP
%% =========================================================================
addpath(spm_path);
addpath(fieldtrip_path);
ft_defaults;

if ~exist(save_dir,'dir'), mkdir(save_dir); end

%% =========================================================================
%  LOAD GEOMETRY — Z-labelled back sensors (z <= -45 mm)
%% =========================================================================
fprintf('Loading geometry...\n');
geom    = load(geomfile);
sensors = geom.experimental_sensors;

is_back   = sensors.chanpos(:,3) <= -45;
is_Z      = strncmp(sensors.label, 'Z', 1);
keep      = is_back & is_Z;

radial_idx        = find(keep);
sensor_pos_radial = sensors.chanpos(radial_idx, :);
radial_labels     = sensors.label(radial_idx);

fprintf('  Back Z-channel sensors: %d\n', numel(radial_idx));

%% =========================================================================
%  LOAD BEM LEADFIELD
%% =========================================================================
fprintf('Loading BEM leadfield...\n');
lf_bem_data = load(lf_bem_path);
lf_bem_raw  = lf_bem_data.leadfield_cord;

pos_bem    = lf_bem_raw.pos * 1000;   % metres -> mm
n_src_bem  = numel(lf_bem_raw.leadfield);
cord_pos_bem = pos_bem(:,2);          % y = inferior->superior

% Match Z-channel rows
[~, bem_lf_rows, bem_sens_rows] = intersect(lf_bem_raw.label, radial_labels, 'stable');
sensor_pos_bem = sensor_pos_radial(bem_sens_rows, :);
fprintf('  BEM: %d sources, %d Z-channel rows matched\n', n_src_bem, numel(bem_lf_rows));

%% =========================================================================
%  LOAD BSLAW LEADFIELD
%% =========================================================================
fprintf('Loading BSLaw leadfield...\n');
lf_bs_data = load(lf_bs_path);
lf_bs_raw  = lf_bs_data.leadfield_bs;

pos_bs      = lf_bs_raw.pos;          % already mm
n_src_bs    = numel(lf_bs_raw.leadfield);
cord_pos_bs = pos_bs(:,2);

[~, bs_lf_rows, bs_sens_rows] = intersect(lf_bs_raw.label, radial_labels, 'stable');
sensor_pos_bs = sensor_pos_radial(bs_sens_rows, :);
fprintf('  BSLaw: %d sources, %d Z-channel rows matched\n', n_src_bs, numel(bs_lf_rows));

%% =========================================================================
%  COMPUTE RMS AMPLITUDE AND CORRELATION ALONG SPINE
%  For each source point and each orientation, compute:
%    - RMS across Z-channel sensors (amplitude proxy)
%    - Pearson correlation between BEM and BSLaw sensor patterns
%  BSLaw has more sources — match each BSLaw source to nearest BEM source
%  for correlation, but plot RMS independently on own source axis
%% =========================================================================
fprintf('Computing amplitude and correlation along spine...\n');

n_ori = numel(ori_cols);

% RMS per source per orientation
rms_bem = zeros(n_src_bem, n_ori);
rms_bs  = zeros(n_src_bs,  n_ori);

for s = 1:n_src_bem
    if ~isempty(lf_bem_raw.leadfield{s})
        lf = lf_bem_raw.leadfield{s}(bem_lf_rows, :);
        for oi = 1:n_ori
            rms_bem(s,oi) = rms(lf(:, ori_cols(oi)));
        end
    end
end

for s = 1:n_src_bs
    if ~isempty(lf_bs_raw.leadfield{s})
        lf = lf_bs_raw.leadfield{s}(bs_lf_rows, :);
        for oi = 1:n_ori
            rms_bs(s,oi) = rms(lf(:, ori_cols(oi)));
        end
    end
end

% Correlation — match each BEM source to nearest BSLaw source by y position
% then correlate the full sensor pattern (all Z channels, each orientation)
corr_bem_bs = zeros(n_src_bem, n_ori);

for s = 1:n_src_bem
    if isempty(lf_bem_raw.leadfield{s}), continue; end

    % Find nearest BSLaw source
    [~, nearest_bs] = min(abs(cord_pos_bs - cord_pos_bem(s)));
    if isempty(lf_bs_raw.leadfield{nearest_bs}), continue; end

    lf_b = lf_bem_raw.leadfield{s}(bem_lf_rows, :);
    lf_s = lf_bs_raw.leadfield{nearest_bs}(bs_lf_rows, :);

    for oi = 1:n_ori
        r = corrcoef(lf_b(:, ori_cols(oi)), lf_s(:, ori_cols(oi)));
        corr_bem_bs(s, oi) = r(1,2);
    end
end

fprintf('  Done.\n');

%% =========================================================================
%  TOPOPLOT LEADFIELD MATRICES AT source_idx
%% =========================================================================
lf_bem_mat     = lf_bem_raw.leadfield{source_idx}(bem_lf_rows, :);
lf_bs_mat      = lf_bs_raw.leadfield{source_idx}(bs_lf_rows,   :);

% Shared clim per orientation across both models
clims = zeros(n_ori, 1);
for oi = 1:n_ori
    clims(oi) = max([ max(abs(lf_bem_mat(:, ori_cols(oi)))), ...
                      max(abs(lf_bs_mat(:,  ori_cols(oi)))) ]);
end

%% =========================================================================
%  FIGURE 1 — Topoplots: 2 rows x n_ori cols
%% =========================================================================
figure('Color','w','Units','inches','Position',[1 1 5*n_ori 7]);
tiledlayout(2, n_ori, 'TileSpacing','compact', 'Padding','tight');

row_labels   = {sprintf('BEM (src %d)',   source_idx), ...
                sprintf('BSLaw (src %d)', source_idx)};
lf_mats      = {lf_bem_mat,    lf_bs_mat};
sensor_poses = {sensor_pos_bem, sensor_pos_bs};

for mi = 1:2
    for oi = 1:n_ori
        nexttile((mi-1)*n_ori + oi);
        plot_topoplot_meg(sensor_poses{mi}, lf_mats{mi}(:, ori_cols(oi)), ...
            [-clims(oi), clims(oi)]);
        if mi == 1
            title(ori_names{oi}, 'FontWeight','bold','FontSize',12);
        end
        if oi == 1
            ylabel(row_labels{mi}, 'FontWeight','bold','FontSize',11);
        end
    end
end

sgtitle(sprintf('Leadfield topoplots — back Z-sensors — source %d', source_idx), ...
    'FontSize', 12, 'FontWeight','bold');

exportgraphics(gcf, fullfile(save_dir, sprintf('lf_topoplot_source%d.png', source_idx)), ...
    'Resolution', 600);
saveas(gcf, fullfile(save_dir, sprintf('lf_topoplot_source%d.fig', source_idx)));

%% =========================================================================
%  FIGURE 2 — RMS amplitude along spine
%% =========================================================================
figure('Color','w','Position',[100 100 900 400*n_ori]);
tiledlayout(n_ori, 1, 'TileSpacing','compact', 'Padding','tight');

for oi = 1:n_ori
    nexttile;
    hold on;
    plot(cord_pos_bem, rms_bem(:,oi), 'b-',  'LineWidth', 2, 'DisplayName', 'BEM');
    plot(cord_pos_bs,  rms_bs(:,oi),  'r--', 'LineWidth', 2, 'DisplayName', 'BSLaw');
    xline(cord_pos_bem(source_idx), 'k:', 'LineWidth', 1.2, 'HandleVisibility','off');
    ylabel('RMS (fT/nAm)');
    title(ori_names{oi}, 'FontWeight','normal','FontSize',11);
    legend('Location','best','FontSize',10);
    grid on; box on;
    if oi == n_ori
        xlabel('Position along cord (mm, inferior \rightarrow superior)');
    else
        set(gca,'XTickLabel',[]);
    end
end

sgtitle('RMS amplitude along spine — back Z-sensors', ...
    'FontWeight','normal','FontSize',12);

exportgraphics(gcf, fullfile(save_dir,'lf_rms_along_spine.png'),'Resolution',300);
saveas(gcf, fullfile(save_dir,'lf_rms_along_spine.fig'));

%% =========================================================================
%  FIGURE 3 — BEM vs BSLaw correlation along spine
%% =========================================================================
figure('Color','w','Position',[100 100 900 400*n_ori]);
tiledlayout(n_ori, 1, 'TileSpacing','compact', 'Padding','tight');

for oi = 1:n_ori
    nexttile;
    hold on;
    plot(cord_pos_bem, corr_bem_bs(:,oi), 'k-', 'LineWidth', 2);
    xline(cord_pos_bem(source_idx), 'b:', 'LineWidth', 1.2);
    yline(0, 'k:', 'HandleVisibility','off');
    ylim([-1 1]);
    ylabel('Pearson r');
    title(sprintf('%s — BEM vs BSLaw', ori_names{oi}), ...
        'FontWeight','normal','FontSize',11);
    grid on; box on;
    if oi == n_ori
        xlabel('Position along cord (mm, inferior \rightarrow superior)');
    else
        set(gca,'XTickLabel',[]);
    end
end

sgtitle('Sensor pattern correlation along spine — BEM vs BSLaw', ...
    'FontWeight','normal','FontSize',12);

exportgraphics(gcf, fullfile(save_dir,'lf_correlation_along_spine.png'),'Resolution',300);
saveas(gcf, fullfile(save_dir,'lf_correlation_along_spine.fig'));

fprintf('\nDone. Saved to %s\n', save_dir);

%% =========================================================================
%  LOCAL FUNCTION
%% =========================================================================
function plot_topoplot_meg(sensor_pos, values, clim)

    if numel(values) ~= size(sensor_pos, 1)
        error('Values (%d) does not match sensors (%d)', ...
              numel(values), size(sensor_pos, 1));
    end

    x = sensor_pos(:, 1);
    y = sensor_pos(:, 2);

    xi = linspace(min(x), max(x), 150);
    yi = linspace(min(y), max(y), 150);
    [Xi, Yi] = meshgrid(xi, yi);
    Zi = griddata(x, y, values, Xi, Yi, 'v4');

    contourf(Xi, Yi, Zi, 25, 'LineStyle','none');
    hold on;

    data_max = max(abs(Zi(:)), [], 'omitnan');
    if data_max > 0
        levels = linspace(-data_max, data_max, 10);
        levels = levels(2:end-1);
        contour(Xi, Yi, Zi, levels, 'LineColor',[0.15 0.15 0.15], 'LineWidth',0.5);
    end

    scatter(x, y, 4, 'k', 'filled', 'MarkerFaceAlpha',0.25, 'MarkerEdgeColor','none');

    n_cols    = 256; half = n_cols/2;
    blue_half = [linspace(0,1,half)', linspace(0,1,half)', ones(half,1)];
    red_half  = [ones(half,1), linspace(1,0,half)', linspace(1,0,half)'];
    colormap(gca, [blue_half; red_half]);

    cb = colorbar; cb.Label.String = 'Leadfield (fT/nAm)';
    caxis(clim);
    axis equal tight off;
    set(gca, 'XTick',[], 'YTick',[], 'FontSize',16);
    hold off;
end