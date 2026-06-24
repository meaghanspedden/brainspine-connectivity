%% figure_BEM_vs_BS_P1.m
%
% Two-panel comparison: BEM vs BS forward model, Participant 1 (OP00212).
% Each panel shows invp_smooth thresholded by significance on the spinal
% cord mesh with vertebral column (bone) only — no torso/lungs/heart.
%
% Shared colour limits across both panels for direct comparison.
% Saved as two identically-sized PNGs for side-by-side figure assembly.

clear; close all; clc;

%% =========================================================================
%  CONFIG
%% =========================================================================
fieldtrip_path = 'C:\Users\mspedden\Documents\fieldtrip';
spm_path       = 'C:\Users\mspedden\Documents\spm';
save_dir       = 'C:\Users\mspedden\Documents\brainspine_savetest';
geomfile       = 'C:\Leadfields meshes\geometries_experimental.mat';
geomfile_brain = 'C:\Leadfields meshes\geometries_cervical_realistic.mat';

addpath(spm_path); spm('defaults','EEG');
addpath(fieldtrip_path); ft_defaults;
ncol = 256;
addpath(fullfile(fieldtrip_path,'external','matplotlib'));

fig_dir = fullfile(save_dir, 'figures');
if ~exist(fig_dir,'dir'), mkdir(fig_dir); end

invpthr = -log10(0.05);

%% =========================================================================
%  LOAD GEOMETRY
%% =========================================================================
fprintf('Loading geometry...\n');
geom_exp     = load(geomfile);
sources_cent = geom_exp.sources_cent;
mesh_wm      = geom_exp.mesh_wm; mesh_wm.unit = 'mm';
mesh_bone    = geom_exp.mesh_bone;

geom_brain = load(geomfile_brain, 'mesh_brain');
mesh_brain = geom_brain.mesh_brain;

%% =========================================================================
%  LOAD P1 RESULTS — BEM AND BS
%% =========================================================================
bem = load(fullfile(save_dir, '\legacy\subResult_subOP00212_true_bem.mat'));
bs  = load(fullfile(save_dir, 'subResult_subOP00212_BS.mat'));

mask_bem        = logical(bem.mask);
invp_bem        = bem.invp_smooth;
invp_bem(~mask_bem) = 0;

mask_bs         = logical(bs.mask);
invp_bs         = bs.invp_smooth;
invp_bs(~mask_bs) = 0;

fprintf('P1 BEM significant sources: %d\n', sum(mask_bem));
fprintf('P1 BS  significant sources: %d\n', sum(mask_bs));

%% =========================================================================
%  SHARED COLOUR LIMITS
%  Lower = invpthr (p=0.05), upper = max across both models
%% =========================================================================
max_bem = max(invp_bem);
max_bs  = max(invp_bs);
% clim_shared = [invpthr, max(max_bem, max_bs)];
% fprintf('Shared colour limits: [%.3f  %.3f]\n', clim_shared);

%% =========================================================================
%  COLORMAP
%% =========================================================================
cmap_hot = [[0.92 0.92 0.92]; flipud(magma(ncol-1))];


spine_int_bem = interp_to_mesh(invp_bem, sources_cent, mesh_wm);
spine_int_bs  = interp_to_mesh(invp_bs,  sources_cent, mesh_wm);

%% =========================================================================
%  SHARED PLOT CONFIG
%% =========================================================================
cfg_plot = [];
cfg_plot.figure       = 'gcf';
cfg_plot.method       = 'surface';
cfg_plot.funparameter = 'coh';
cfg_plot.funcolormap  = cmap_hot;
%cfg_plot.funcolorlim  = clim_shared;
cfg_plot.projmethod   = 'nearest';
cfg_plot.surffile     = mesh_wm;

fig_sz = [100 100 420 620];   % identical size for both panels

%% =========================================================================
%  PANEL 1 — BEM
%% =========================================================================
hfig_bem = figure('Color','w','Position', fig_sz);
ft_sourceplot(cfg_plot, spine_int_bem);
view(-250,-1); camlight; ax = gca; ax.FontSize = 13;
hold on;
ft_plot_mesh(mesh_brain,'facecolor',[0.8 0.3 0.3],'facealpha',0.07,'edgecolor','none');
ft_plot_mesh(mesh_bone, 'facecolor',[0.9 0.85 0.7],'facealpha',0.35,'edgecolor','none');
title('Participant 1 — BEM','Interpreter','none','FontSize',13,'FontWeight','normal');

set(hfig_bem,'Units','Inches');
pos = get(hfig_bem,'Position');
set(hfig_bem,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3) pos(4)]);
print(hfig_bem, fullfile(fig_dir,'figure_P1_BEM.png'), '-dpng', '-r300');
savefig(hfig_bem, fullfile(fig_dir,'figure_P1_BEM.fig'));
fprintf('BEM panel saved.\n');

%% =========================================================================
%  PANEL 2 — BS
%% =========================================================================
hfig_bs = figure('Color','w','Position', fig_sz);
ft_sourceplot(cfg_plot, spine_int_bs);
view(-250,-1); camlight; ax = gca; ax.FontSize = 13;
hold on;
ft_plot_mesh(mesh_brain,'facecolor',[0.8 0.3 0.3],'facealpha',0.07,'edgecolor','none');
ft_plot_mesh(mesh_bone, 'facecolor',[0.9 0.85 0.7],'facealpha',0.35,'edgecolor','none');
title('Participant 1 — BS','Interpreter','none','FontSize',13,'FontWeight','normal');

set(hfig_bs,'Units','Inches');
pos = get(hfig_bs,'Position');
set(hfig_bs,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3) pos(4)]);
print(hfig_bs, fullfile(fig_dir,'figure_P1_BS.png'), '-dpng', '-r300');
savefig(hfig_bs, fullfile(fig_dir,'figure_P1_BS.fig'));
fprintf('BS panel saved.\n');

fprintf('\n=== DONE ===\n');
fprintf('Colour limits: [%.3f  %.3f]  (shared across both panels)\n', clim_shared);

%% =========================================================================
%  INTERPOLATE BOTH ONTO MESH
%% =========================================================================
function spine_int = interp_to_mesh(invp_display, sources_cent, mesh_wm)
    source_p         = [];
    source_p.pos     = sources_cent.pos;
    source_p.inside  = sources_cent.inside;
    source_p.avg.coh = invp_display;
    cfg_interp             = [];
    cfg_interp.parameter   = 'coh';
    spine_int = ft_sourceinterpolate(cfg_interp, source_p, mesh_wm);
end
