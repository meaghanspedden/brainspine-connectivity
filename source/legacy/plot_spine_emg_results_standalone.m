%% plot_spine_results_standalone.m
% Loads saved spine-EMG DICS results and produces:
%   - Summary of significant participants
%   - P1: mesh-only and compound figure (invp_smooth)
%   - Group: mesh-only and compound prevalence figure
%   - Overlap: P1 coherence vs group prevalence line plot

clear all; close all; clc;

%% =========================================================================
%  CONFIG
%% =========================================================================
fieldtrip_path = 'C:\Users\mspedden\Documents\fieldtrip';
save_dir       = 'C:\Users\mspedden\Documents\brainspine_savetest';
geomfile       = 'C:\Leadfields meshes\geometries_experimental.mat';
geomfile_brain = 'C:\Leadfields meshes\geometries_cervical_realistic.mat';
fig_dir        = fullfile(save_dir, 'figures');
out_suffix     = '_BS';
threshold      = 0.2;

subs_spine = {'OP00212','OP00213','OP00215','OP00219', ...
              'OP00220','OP00221','OP00224','OP00225','OP00226'};

%% =========================================================================
%  SETUP
%% =========================================================================
addpath(fieldtrip_path);
ft_defaults;
addpath(fullfile(fieldtrip_path,'external','matplotlib'));
if ~exist(fig_dir,'dir'), mkdir(fig_dir); end

ncol     = 256;
cmap_hot = [[0.92 0.92 0.92]; flipud(magma(ncol-1))];
invpthr  = -log10(0.05);

%% =========================================================================
%  LOAD GEOMETRY
%% =========================================================================
fprintf('Loading geometry...\n');
geom_exp      = load(geomfile);
sources_cent  = geom_exp.sources_cent;
mesh_torso    = geom_exp.mesh_torso;
mesh_wm       = geom_exp.mesh_wm;
mesh_bone     = geom_exp.mesh_bone;
mesh_wm.unit  = 'mm';
nsourcepoints = size(sources_cent.pos, 1);

geom_brain = load(geomfile_brain, 'mesh_brain');
mesh_brain = geom_brain.mesh_brain;
mesh_cut   = clip_torso(mesh_torso);

%% =========================================================================
%  LOAD GROUP RESULTS
%% =========================================================================
fprintf('Loading group results...\n');
group_file = fullfile(save_dir, ['groupRes_spine_DICS_' out_suffix '.mat']);
assert(exist(group_file,'file')==2, 'Group results file not found: %s', group_file);
loaded      = load(group_file, 'subjResults');
subjResults = loaded.subjResults;
nSubjects   = numel(subjResults);

%% =========================================================================
%  SUMMARY
%% =========================================================================
fprintf('\n── Significant spine-EMG coherence ──\n');
sig_pos = false(nSubjects, 1);
for ss = 1:nSubjects
    sig_pos(ss) = any(subjResults(ss).coh_diff > subjResults(ss).thr95);
    fprintf('  %s: %s\n', subs_spine{ss}, mat2str(sig_pos(ss)));
end
fprintf('  Total: %d / %d\n\n', sum(sig_pos), nSubjects);

%% =========================================================================
%  GROUP PREVALENCE MAP
%% =========================================================================
all_masks = zeros(nsourcepoints, nSubjects);
for s = 1:nSubjects
    all_masks(:,s) = double(subjResults(s).coh_diff > subjResults(s).thr95);
end
prevalence_loc = mean(all_masks, 2);
gp_masked      = prevalence_loc;
gp_masked(gp_masked < threshold) = 0;

group_ft        = [];
group_ft.pos    = sources_cent.pos;
group_ft.inside = sources_cent.inside;
group_ft.pow    = gp_masked;

cfg_gi    = []; cfg_gi.parameter = 'pow'; cfg_gi.interpmethod = 'nearest';
group_int = ft_sourceinterpolate(cfg_gi, group_ft, mesh_wm);
clim_grp  = [threshold max(group_int.pow)];

prev_vals = gp_masked; prev_vals(prev_vals < threshold) = -inf;
has_grp   = any(isfinite(prev_vals) & prev_vals > -inf);
if has_grp
    peaks_grp = sources_cent.pos(prev_vals >= max(prev_vals(isfinite(prev_vals)))*0.99, :);
    fprintf('  Group peak cranio-caudal position: %.1f mm\n', mean(peaks_grp(:,2)));
end

cfg_grp               = [];
cfg_grp.method        = 'surface';
cfg_grp.funparameter  = 'pow';
cfg_grp.funcolorlim   = clim_grp;
cfg_grp.funcolormap   = cmap_hot;
cfg_grp.projmethod    = 'nearest';
cfg_grp.surffile      = mesh_wm;
cfg_grp.opacitylim    = clim_grp;
cfg_grp.opacitymap    = 'rampup';

%% ── Group mesh only ──
hfig_grp_mesh        = figure('Color','w','Position',[100 100 400 650]);
cfg_grp.figure       = 'gcf';
ft_sourceplot(cfg_grp, group_int);
colorbar off;
view(-250,-1); camlight; ax = gca; ax.FontSize = 14;
% title('Group prevalence — spine-EMG (BS)', 'Interpreter','none', 'FontSize', 13);
% if has_grp
%     hold on;
%     scatter3(peaks_grp(:,1), peaks_grp(:,2), peaks_grp(:,3)+10, ...
%         50, 'o', 'filled', 'MarkerFaceColor','k', 'MarkerEdgeColor','k');
%     uistack(findobj(gca,'Type','Scatter'), 'top');
% end
print(hfig_grp_mesh, fullfile(fig_dir, ['group_spine_meshonly' out_suffix '.png']), '-dpng', '-r300');

%% ── Group compound ──
hfig_grp_comp  = figure('Color','w');
cfg_grp.figure = 'gcf';
ft_sourceplot(cfg_grp, group_int);
view(-250,-1); camlight; ax = gca; ax.FontSize = 14;
hold on;
ft_plot_mesh(mesh_brain, 'facecolor',[0.8 0.3 0.3], 'facealpha',0.07, 'edgecolor','none');
ft_plot_mesh(mesh_cut,   'facecolor',[0.3 0.3 0.9], 'facealpha',0.10, 'edgecolor','none');
ft_plot_mesh(mesh_bone,  'facecolor',[0.9 0.85 0.7],'facealpha',0.30, 'edgecolor','none');
% if has_grp
%     scatter3(peaks_grp(:,1), peaks_grp(:,2), peaks_grp(:,3)+10, ...
%         50, 'o', 'filled', 'MarkerFaceColor','k', 'MarkerEdgeColor','k');
%     uistack(findobj(gca,'Type','Scatter'), 'top');
% end
%print(hfig_grp_comp, fullfile(fig_dir, ['group_spine_compound' out_suffix '.png']), '-dpng', '-r300');

%% =========================================================================
%  P1 FIGURES
%% =========================================================================
fprintf('Plotting P1 figures...\n');
p1          = subjResults(1);
invp_smooth = p1.invp_smooth;
mask_p1     = logical(p1.sig_mask);

source_p         = [];
source_p.pos     = sources_cent.pos;
source_p.inside  = true(nsourcepoints, 1);
source_p.avg.coh = invp_smooth;
if isfield(sources_cent, 'dim')
    source_p.dim = sources_cent.dim;
end

cfg_interp = []; cfg_interp.parameter = 'coh';
spine_int  = ft_sourceinterpolate(cfg_interp, source_p, mesh_wm);

invp_max = max(invp_smooth);
if invp_max <= invpthr
    clim_p1 = [invpthr invpthr + 0.5];
else
    clim_p1 = [invpthr invp_max];
end

sig_invp = invp_smooth; sig_invp(~mask_p1) = -inf;
has_sig  = any(isfinite(sig_invp) & sig_invp > -inf);
if has_sig
    peaks_p1 = sources_cent.pos(sig_invp >= max(sig_invp(isfinite(sig_invp)))*0.99, :);
end

cfg_plot              = [];
cfg_plot.method       = 'surface';
cfg_plot.funparameter = 'coh';
cfg_plot.funcolormap  = cmap_hot;
cfg_plot.funcolorlim  = clim_p1;
cfg_plot.projmethod   = 'nearest';
cfg_plot.surffile     = mesh_wm;

%% ── P1 mesh only ──
hfig_p1_mesh   = figure('Color','w');
cfg_plot.figure = 'gcf';
ft_sourceplot(cfg_plot, spine_int);
colorbar off;
view(-250,-1); camlight; ax = gca; ax.FontSize = 14;
%title('Participant 1 — spine-EMG coherence (BS)', 'Interpreter','none', 'FontSize', 13);
% if has_sig
%     hold on;
%     scatter3(peaks_p1(:,1), peaks_p1(:,2), peaks_p1(:,3)+10, ...
%         50, 'o', 'filled', 'MarkerFaceColor','k', 'MarkerEdgeColor','k');
%     uistack(findobj(gca,'Type','Scatter'), 'top');
% end
%print(hfig_p1_mesh, fullfile(fig_dir, ['P1_spine_meshonly' out_suffix '.png']), '-dpng', '-r300');

%% ── P1 compound ──
hfig_p1_comp   = figure('Color','w');
cfg_plot.figure = 'gcf';
ft_sourceplot(cfg_plot, spine_int);
view(-250,-1); camlight; ax = gca; ax.FontSize = 14;
hold on;
ft_plot_mesh(mesh_brain, 'facecolor',[0.8 0.3 0.3], 'facealpha',0.07, 'edgecolor','none');
ft_plot_mesh(mesh_cut,   'facecolor',[0.3 0.3 0.9], 'facealpha',0.10, 'edgecolor','none');
ft_plot_mesh(mesh_bone,  'facecolor',[0.9 0.85 0.7],'facealpha',0.30, 'edgecolor','none');
% if has_sig
%     scatter3(peaks_p1(:,1), peaks_p1(:,2), peaks_p1(:,3)+10, ...
%         50, 'o', 'filled', 'MarkerFaceColor','k', 'MarkerEdgeColor','k');
%     uistack(findobj(gca,'Type','Scatter'), 'top');
% end
% title('Participant 1 — spine-EMG coherence (BS)', 'Interpreter','none');
% print(hfig_p1_comp, fullfile(fig_dir, ['P1_spine_compound' out_suffix '.png']), '-dpng', '-r300');

%% =========================================================================
%  OVERLAP ON MESH
%% =========================================================================
%% =========================================================================
%  OVERLAP ON MESH — overlap only, compound view
%% =========================================================================
% Binary map: 1 = overlap only (P1 sig AND group prevalence >= threshold)
overlap_map = zeros(nsourcepoints, 1);
overlap_map(logical(p1.sig_mask) & prevalence_loc >= 0.3) = 1;

% Interpolate onto mesh
source_ov         = [];
source_ov.pos     = sources_cent.pos;
source_ov.inside  = true(nsourcepoints, 1);
source_ov.avg.coh = overlap_map;
cfg_oi            = []; cfg_oi.parameter = 'coh';
ov_int            = ft_sourceinterpolate(cfg_oi, source_ov, mesh_wm);

% Two-colour map: grey background, teal for overlap
cmap_ov = [0.92 0.92 0.92;    % 0: background
           0.08 0.70 0.68];   % 1: overlap (teal)

hfig_ov_mesh = figure('Color','w','Position',[100 100 500 650]);
cfg_ov              = [];
cfg_ov.figure       = 'gcf';
cfg_ov.method       = 'surface';
cfg_ov.funparameter = 'coh';
cfg_ov.funcolormap  = cmap_ov;
cfg_ov.funcolorlim  = [0 1];
cfg_ov.projmethod   = 'nearest';
cfg_ov.surffile     = mesh_wm;
ft_sourceplot(cfg_ov, ov_int);
colorbar off;
view(-250,-1); camlight; ax = gca; ax.FontSize = 14;
hold on;

% Compound overlays
ft_plot_mesh(mesh_brain, 'facecolor',[0.8 0.3 0.3], 'facealpha',0.07, 'edgecolor','none');
ft_plot_mesh(mesh_cut,   'facecolor',[0.3 0.3 0.9], 'facealpha',0.10, 'edgecolor','none');
ft_plot_mesh(mesh_bone,  'facecolor',[0.9 0.85 0.7],'facealpha',0.30, 'edgecolor','none');

title(sprintf('P1 \\cap group (prev\\geq%.1f)', 0.3), ...
    'FontSize', 13, 'FontWeight','normal');



%% =========================================================================
%  HELPER
%% =========================================================================
function mesh_cut = clip_torso(mesh_torso)
    y         = mesh_torso.vertices(:,2);
    keep_vert = y > -200;
    new_idx   = zeros(size(keep_vert));
    new_idx(keep_vert) = 1:sum(keep_vert);
    faces_keep        = all(keep_vert(mesh_torso.faces), 2);
    mesh_cut.vertices = mesh_torso.vertices(keep_vert,:);
    mesh_cut.faces    = new_idx(mesh_torso.faces(faces_keep,:));
    mesh_cut.unit     = mesh_torso.unit;
end