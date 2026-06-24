%% figure_BEM_vs_BS_P1.m
%
% Two individual figures: BEM vs BS, Participant 1 (OP00212).
% invp_smooth on spinal cord mesh + vertebral column only.
% Shared colour limits. Colorbar on each figure.

clear; close all; clc;

%% =========================================================================
%  CONFIG
%% =========================================================================
fieldtrip_path = 'C:\Users\mspedden\Documents\fieldtrip';
spm_path       = 'C:\Users\mspedden\Documents\spm';
save_dir       = 'C:\Users\mspedden\Documents\brainspine_savetest';
geomfile       = 'C:\Leadfields meshes\geometries_experimental.mat';

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

%% =========================================================================
%  LOAD P1 RESULTS
%% =========================================================================
fprintf('Loading P1 BEM results...\n');
bem = load(fullfile(save_dir, 'legacy', 'subResult_subOP00212_true_bem.mat'));

fprintf('Loading P1 BS results...\n');
bs  = load(fullfile(save_dir, 'subResult_subOP00212_BS.mat'));

mask_bem = logical(bem.mask);
invp_bem = bem.invp_smooth;

mask_bs  = logical(bs.mask);
invp_bs  = bs.invp_smooth;

fprintf('P1 BEM significant sources: %d / %d\n', sum(mask_bem), numel(mask_bem));
fprintf('P1 BS  significant sources: %d / %d\n', sum(mask_bs),  numel(mask_bs));

%% =========================================================================
%  SHARED COLOUR LIMITS
%% =========================================================================
clim_shared = [invpthr, max([max(invp_bem), max(invp_bs), invpthr + 0.5])];
fprintf('Shared colour limits: [%.3f  %.3f]\n', clim_shared);

%% =========================================================================
%  COLORMAP
%% =========================================================================
cmap_hot = [[0.92 0.92 0.92]; flipud(magma(ncol-1))];

%% =========================================================================
%  INTERPOLATE ONTO MESH
%% =========================================================================
source_bem         = [];
source_bem.pos     = sources_cent.pos;
source_bem.inside  = logical(sources_cent.inside);
source_bem.unit    = 'mm';
source_bem.avg.coh = invp_bem;

cfg_interp           = [];
cfg_interp.parameter = 'coh';
spine_int_bem = ft_sourceinterpolate(cfg_interp, source_bem, mesh_wm);

source_bs         = [];
source_bs.pos     = sources_cent.pos;
source_bs.inside  = logical(sources_cent.inside);
source_bs.unit    = 'mm';
source_bs.avg.coh = invp_bs;

spine_int_bs = ft_sourceinterpolate(cfg_interp, source_bs, mesh_wm);

%% =========================================================================
%  FIGURE 1 — BEM
%% =========================================================================
hfig_bem = figure('Color','w','Position',[100 100 480 640]);
cfg_plot = [];
cfg_plot.figure       = 'gcf';
cfg_plot.method       = 'surface';
cfg_plot.funparameter = 'coh';
cfg_plot.funcolormap  = cmap_hot;
cfg_plot.funcolorlim  = clim_shared;
cfg_plot.projmethod   = 'nearest';
cfg_plot.surffile     = mesh_wm;
ft_sourceplot(cfg_plot, spine_int_bem);
delete(findobj(hfig_bem, 'Type', 'colorbar'));
view(-250,-1); camlight; ax = gca; ax.FontSize = 13; hold on;
ft_plot_mesh(mesh_bone,'facecolor',[0.9 0.85 0.7],'facealpha',0.35,'edgecolor','none');
cb = colorbar(ax,'Location','southoutside');
cb.Label.String   = 'Contraction > Rest coherence (-log_{10}(p))';
cb.Label.FontSize = 11;
clim(ax, clim_shared);
title('Participant 1 — BEM','FontSize',13,'FontWeight','normal','Interpreter','none');

set(hfig_bem,'Units','Inches');
pos = get(hfig_bem,'Position');
set(hfig_bem,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3) pos(4)]);
print(hfig_bem, fullfile(fig_dir,'figure_P1_BEM.png'),'-dpng','-r300');
savefig(hfig_bem, fullfile(fig_dir,'figure_P1_BEM.fig'));
fprintf('BEM figure saved.\n');

%% =========================================================================
%  FIGURE 2 — BS
%% =========================================================================
hfig_bs = figure('Color','w','Position',[600 100 480 640]);
ft_sourceplot(cfg_plot, spine_int_bs);
delete(findobj(hfig_bs, 'Type', 'colorbar'));
view(-250,-1); camlight; ax = gca; ax.FontSize = 13; hold on;
ft_plot_mesh(mesh_bone,'facecolor',[0.9 0.85 0.7],'facealpha',0.35,'edgecolor','none');
cb = colorbar(ax,'Location','southoutside');
cb.Label.String   = 'Contraction > Rest coherence (-log_{10}(p))';
cb.Label.FontSize = 11;
clim(ax, clim_shared);
title('Participant 1 — BS','FontSize',13,'FontWeight','normal','Interpreter','none');

set(hfig_bs,'Units','Inches');
pos = get(hfig_bs,'Position');
set(hfig_bs,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3) pos(4)]);
print(hfig_bs, fullfile(fig_dir,'figure_P1_BS.png'),'-dpng','-r300');
savefig(hfig_bs, fullfile(fig_dir,'figure_P1_BS.fig'));
fprintf('BS figure saved.\n');

fprintf('\n=== DONE ===\n');
fprintf('BEM sig: %d  |  BS sig: %d  |  shared clim=[%.2f %.2f]\n', ...
    sum(mask_bem), sum(mask_bs), clim_shared);