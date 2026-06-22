%% plot_p1_spine_emg.m
% Load saved results for OP00212 and reproduce spine-EMG coherence figures.

clear all; close all; clc;

%% =========================================================================
%  USER CONFIG
%% =========================================================================
fieldtrip_path = 'C:\Users\mspedden\Documents\fieldtrip';
spm_path       = 'C:\Users\mspedden\Documents\spm';
save_dir       = 'C:\Users\mspedden\Documents\brainspine_savetest';
geomfile       = 'C:\Leadfields meshes\geometries_experimental.mat';
geomfile_brain = 'C:\Leadfields meshes\geometries_cervical_realistic.mat';

out_suffix  = '_true_bem';
sub         = 'OP00212';
fwhm_mm     = 20;

roi_sensory = 14:20;
roi_motor   = 21:27;
roi_all     = 14:27;
z_offset    = 10;

%% =========================================================================
%  SETUP
%% =========================================================================
addpath(spm_path);
addpath(fieldtrip_path);
ft_defaults;

fig_dir = fullfile(save_dir, 'figures');
if ~exist(fig_dir,'dir'), mkdir(fig_dir); end

%% =========================================================================
%  LOAD GEOMETRY
%% =========================================================================
fprintf('Loading geometry...\n');
geom_exp     = load(geomfile);
sources_cent = geom_exp.sources_cent;
mesh_torso   = geom_exp.mesh_torso;
mesh_wm      = geom_exp.mesh_wm;
mesh_bone    = geom_exp.mesh_bone;
mesh_lungs   = geom_exp.mesh_lungs;
mesh_heart   = geom_exp.mesh_heart;
mesh_wm.unit = 'mm';

geom_brain = load(geomfile_brain, 'mesh_brain');
mesh_brain = geom_brain.mesh_brain;

cord_pos      = sources_cent.pos(:,2);
nsourcepoints = size(sources_cent.pos, 1);
mesh_cut      = clip_torso(mesh_torso);

%% =========================================================================
%  COLOURMAP
%% =========================================================================
ncol     = 256;
addpath(fullfile(fieldtrip_path,'external','matplotlib'));
cmap_hot = [[0.92 0.92 0.92]; flipud(magma(ncol-1))];

%% =========================================================================
%  LOAD RESULTS
%% =========================================================================
fprintf('Loading results for %s...\n', sub);
res = load(fullfile(save_dir, sprintf('subResult_sub%s%s.mat', sub, out_suffix)));

coh_diff     = res.coh_diff;
cohDiff_perm = res.cohDiff_perm;
thr95        = res.thr95;
mask         = res.mask;
invp_smooth  = res.invp_smooth;
nPerm        = size(cohDiff_perm, 2);

fprintf('  thr95=%.6f  sig sources=%d/%d  nPerm=%d\n', ...
    thr95, sum(mask), nsourcepoints, nPerm);

%% =========================================================================
%  INTERPOLATE ONTO SPINE MESH
%% =========================================================================
source_p         = [];
source_p.pos     = sources_cent.pos;
source_p.inside  = true(nsourcepoints, 1);
source_p.avg.coh = invp_smooth;
if isfield(sources_cent, 'dim')
    source_p.dim = sources_cent.dim;
end

cfg_interp = []; cfg_interp.parameter = 'coh';
spine_int  = ft_sourceinterpolate(cfg_interp, source_p, mesh_wm);

invpthr  = -log10(0.05);
invp_max = max(invp_smooth);
if invp_max <= invpthr
    clim_spine = [invpthr invpthr + 0.5];
else
    clim_spine = [invpthr invp_max];
end

%% =========================================================================
%  FIGURE 1 — Spine mesh with torso
%% =========================================================================
fprintf('Figure 1: spine mesh with torso...\n');
hfig_spine = figure;
cfg_plot = []; cfg_plot.figure = 'gcf'; cfg_plot.method = 'surface';
cfg_plot.funparameter = 'coh'; cfg_plot.funcolormap = cmap_hot;
cfg_plot.funcolorlim  = clim_spine; cfg_plot.projmethod = 'nearest';
cfg_plot.surffile = mesh_wm;
ft_sourceplot(cfg_plot, spine_int);
view(-250,-1); camlight; ax = gca; ax.FontSize = 14;
hold on;
ft_plot_mesh(mesh_brain,'facecolor',[0.8 0.3 0.3],'facealpha',0.07,'edgecolor','none');
ft_plot_mesh(mesh_cut,  'facecolor',[0.3 0.3 0.9],'facealpha',0.1, 'edgecolor','none');
ft_plot_mesh(mesh_bone, 'facecolor',[0.9 0.85 0.7],'facealpha',0.3, 'edgecolor','none');
ft_plot_mesh(mesh_lungs,'facecolor',[0.8 0.3 0.3],'facealpha',0.1, 'edgecolor','none');
ft_plot_mesh(mesh_heart,'facecolor',[0.8 0.3 0.3],'facealpha',0.1, 'edgecolor','none');
add_roi_markers(sources_cent, roi_sensory, roi_motor, z_offset);
title(sprintf('%s — spine-EMG coherence (true BEM)', sub),'Interpreter','none');
savefig(hfig_spine, fullfile(fig_dir, ...
    sprintf('step2_sub%s_spineEMG_coherence%s.fig', sub, out_suffix)));
%close(hfig_spine);

%% =========================================================================
%  FIGURE 2 — Mesh only, no colorbar, peak marker + ROI
%% =========================================================================
fprintf('Figure 2: mesh only...\n');
hfig_mesh = figure('Color','w');
cfg_m = []; cfg_m.figure = 'gcf'; cfg_m.method = 'surface';
cfg_m.funparameter = 'coh'; cfg_m.funcolormap = cmap_hot;
cfg_m.funcolorlim  = clim_spine; cfg_m.projmethod = 'nearest';
cfg_m.surffile = mesh_wm;
ft_sourceplot(cfg_m, spine_int);
colorbar off;
view(-250,-1); camlight; ax = gca; ax.FontSize = 14;
title('Participant 1 — spine-EMG coherence (true BEM)', ...
    'Interpreter','none','FontSize',13);
hold on;

% Peak marker
sig_invp = invp_smooth; sig_invp(~mask) = -inf;
if any(isfinite(sig_invp)) && any(sig_invp > -inf)
    peaks_p1 = sources_cent.pos(sig_invp >= ...
        max(sig_invp(isfinite(sig_invp)))*0.99, :);
    scatter3(peaks_p1(:,1), peaks_p1(:,2), peaks_p1(:,3)+z_offset, ...
        200, '*', 'MarkerEdgeColor',[1 1 0], 'LineWidth', 2);
    scatter_obj = findobj(gca,'Type','Scatter');
    uistack(scatter_obj,'top');
end
add_roi_markers(sources_cent, roi_sensory, roi_motor, z_offset);
savefig(hfig_mesh, fullfile(fig_dir, ...
    sprintf('step2_sub%s_spineEMG_meshonly%s.fig', sub, out_suffix)));
%close(hfig_mesh);

%% =========================================================================
%  FIGURE 3 — Coherence diff line plot with ROI shading
%% =========================================================================
fprintf('Figure 3: coherence line plot...\n');
hfig_line = figure('Color','w','Position',[100 100 750 450]);
hold on;

yl_pad = [min(coh_diff)*1.1, max(coh_diff)*1.1];
fill([cord_pos(roi_all(1))   cord_pos(roi_all(end)) ...
      cord_pos(roi_all(end)) cord_pos(roi_all(1))], ...
     [yl_pad(1) yl_pad(1) yl_pad(2) yl_pad(2)], ...
     [0.85 0.85 0.85], 'EdgeColor','none', 'FaceAlpha',0.3, ...
     'DisplayName','ROI (C6-T1)');

plot(cord_pos, coh_diff, 'b-', 'LineWidth', 2, 'DisplayName','Coh diff (stat-rest)');
yline(thr95, 'r--', 'LineWidth', 1.5, ...
    'DisplayName', sprintf('Threshold (%.2e)', thr95));
if any(mask)
    scatter(cord_pos(mask), coh_diff(mask), 60, 'r', 'filled', ...
        'DisplayName','Significant');
end
yline(0,'k:','HandleVisibility','off');

ylim(yl_pad);
xlim([min(cord_pos) max(cord_pos)]);
xlabel('Position along cord (mm, inferior \rightarrow superior)','FontSize',13);
ylabel('Coherence difference (stat-rest)','FontSize',13);
title(sprintf('%s — spine-EMG coherence diff (true BEM, smoothed %d mm)', ...
    sub, fwhm_mm),'Interpreter','none');
legend('Location','northwest','FontSize',10);
grid on; box on; set(gca,'FontSize',13);
savefig(hfig_line, fullfile(fig_dir, ...
    sprintf('step2_sub%s_cohdiff_line%s.fig', sub, out_suffix)));
%close(hfig_line);

%% =========================================================================
%  FIGURE 4 — Null maxima histogram
%% =========================================================================
fprintf('Figure 4: null maxima histogram...\n');
[~, maxIdx_perm] = max(cohDiff_perm, [], 1);
[~, obsMaxIdx]   = max(coh_diff);

hfig_null = figure('Color','w','Position',[100 100 600 450]);
hold on;
histogram(cord_pos(maxIdx_perm), 44, ...
    'FaceColor',[0.75 0.75 0.75],'EdgeColor','k','LineWidth',0.8);
xline(cord_pos(obsMaxIdx), '-', 'Color',[0.2 0 0], 'LineWidth', 2);
xlabel('Cranio-caudal position (mm)','FontSize',14);
ylabel('Count','FontSize',14);
legend({'Null maxima','Observed maximum'},'Location','best','FontSize',14,'Box','off');
set(gca,'FontSize',14,'LineWidth',1.2,'TickDir','out'); box off;
title(sprintf('%s — null maxima (true BEM)', sub),'Interpreter','none');
savefig(hfig_null, fullfile(fig_dir, ...
    sprintf('step2_sub%s_null_maxima%s.fig', sub, out_suffix)));
%close(hfig_null);

%% =========================================================================
%  FIGURE 5 — Null distribution of max coh diff
%% =========================================================================
fprintf('Figure 5: null distribution...\n');
maxPerm     = max(cohDiff_perm, [], 1);
maxPerm_pos = max(maxPerm, 0);

nBins     = 30;
bin_edges = linspace(0, max(maxPerm_pos)*1.05, nBins+1);
bin_ctrs  = (bin_edges(1:end-1) + bin_edges(2:end)) / 2;

hfig_nulldist = figure('Color','w','Position',[100 100 500 430]);
histogram(maxPerm_pos, bin_edges, ...
    'FaceColor',[0.75 0.75 0.75],'EdgeColor','k','LineWidth',0.8);
hold on;

[~, thr_bin] = min(abs(bin_ctrs - thr95));
[~, obs_bin] = min(abs(bin_ctrs - max(coh_diff)));
xline(bin_ctrs(thr_bin), '--k', 'LineWidth', 2,     'DisplayName','Threshold (p<0.05)');
xline(bin_ctrs(obs_bin), '-',   'Color',[0.55 0 0], 'LineWidth', 2.5, ...
    'DisplayName','Observed maximum');

scale_fac = 1e4;
tick_step = max(1, floor(nBins/8));
tick_idx  = 1:tick_step:nBins;
xticks(bin_ctrs(tick_idx));
xticklabels(arrayfun(@(v) sprintf('%.2f', v*scale_fac), bin_ctrs(tick_idx), ...
    'UniformOutput', false));
xtickangle(35);

xlabel(sprintf('Max coherence diff (stat-rest) \\times 10^{-%d}', ...
    round(log10(1/scale_fac))),'FontSize',14);
ylabel('Count','FontSize',14);
legend('Location','best','FontSize',11,'Box','off');
set(gca,'FontSize',13,'LineWidth',1.1,'TickDir','out'); box off;
title(sprintf('%s — null distribution (true BEM, smoothed %d mm)', sub, fwhm_mm), ...
    'Interpreter','none','FontSize',13);
savefig(hfig_nulldist, fullfile(fig_dir, ...
    sprintf('step2_sub%s_null_distribution%s.fig', sub, out_suffix)));
%close(hfig_nulldist);

fprintf('\nDone. All figures saved to %s\n', fig_dir);

%% =========================================================================
%  LOCAL FUNCTIONS
%% =========================================================================
function add_roi_markers(sources_cent, roi_sensory, roi_motor, z_offset)
    for s = roi_sensory
        plot3(sources_cent.pos(s,1), sources_cent.pos(s,2), ...
            sources_cent.pos(s,3) + z_offset, ...
            'o', 'MarkerFaceColor',[0.2 0.4 0.8], 'MarkerEdgeColor','w', ...
            'MarkerSize',10, 'LineWidth',1);
    end
    for s = roi_motor
        plot3(sources_cent.pos(s,1), sources_cent.pos(s,2), ...
            sources_cent.pos(s,3) + z_offset, ...
            'o', 'MarkerFaceColor',[0.8 0.2 0.2], 'MarkerEdgeColor','w', ...
            'MarkerSize',10, 'LineWidth',1);
    end
    plot3(nan,nan,nan,'o','MarkerFaceColor',[0.2 0.4 0.8],'MarkerEdgeColor','w', ...
        'MarkerSize',10,'DisplayName','Sensory ROI (C6-C8)');
    plot3(nan,nan,nan,'o','MarkerFaceColor',[0.8 0.2 0.2],'MarkerEdgeColor','w', ...
        'MarkerSize',10,'DisplayName','Motor ROI (C8-T1)');
    legend('Location','best','FontSize',10);
end

function mesh_cut = clip_torso(mesh_torso)
    y = mesh_torso.vertices(:,2);
    keep_vert = y > -200;
    new_idx   = zeros(size(keep_vert));
    new_idx(keep_vert) = 1:sum(keep_vert);
    faces_keep        = all(keep_vert(mesh_torso.faces), 2);
    mesh_cut.vertices = mesh_torso.vertices(keep_vert,:);
    mesh_cut.faces    = new_idx(mesh_torso.faces(faces_keep,:));
    mesh_cut.unit     = mesh_torso.unit;
end