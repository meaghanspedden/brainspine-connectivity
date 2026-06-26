%% RUN_FIGURES_ONLY.m
% Regenerates all pipeline figures from pre-computed results.
%
% Requires the following saved files in save_dir:
%   groupRes_brain_DICS_brain_pct10.mat
%   groupRes_spine_DICS__BS.mat
%   cluster_spineEMG_pos_BS.mat
%   brainspine_boxplot_data_bslaw_prevalence.mat
%   geometries_experimental_withbrain.mat  (geomfile)
%   T.mat
%
% EDIT THE USER CONFIG SECTION BELOW, THEN RUN.
% =========================================================================

clear all
close all
clc

%% =========================================================================
%  USER CONFIG — edit this section only
%  =========================================================================

% --- Toolbox / path setup -------------------------------------------------
cfg.fieldtrip_path  = 'C:\Users\mspedden\Documents\fieldtrip';   % update if different
cfg.spm_path        = 'C:\Users\mspedden\Documents\spm';          % update if different

% --- Data & geometry paths ------------------------------------------------
cfg.geomfile = 'C:\Leadfields meshes\geometries_experimental_withbrain.mat';
cfg.T_mat_path      = 'C:\Leadfields meshes\T.mat';
cfg.save_dir        = 'C:\forGeneratingFigures';

% --- Smoothing config 
cfg.doSmooth             = 1;
cfg.spine_smooth_fwhm_mm = 20;
cfg.brain_smooth_fwhm_mm = 8;


% --- Steps to plot --------------------------------------------------------
plot_step1 = 1;   % Brain-EMG coherence figures
plot_step2 = 1;   % Spinal cord-EMG coherence figures
plot_step4 = 1;   % Brain-SpineVE coherence + M1 overlap figures

% --- Figure saving --------------------------------------------------------
cfg.saveFigs = 1;

%% =========================================================================
%  SETUP
%% =========================================================================

addpath(cfg.spm_path)
spm('defaults','EEG')
addpath(cfg.fieldtrip_path)
ft_defaults

fig_dir = fullfile(cfg.save_dir, 'figures');
if cfg.saveFigs && ~exist(fig_dir,'dir'), mkdir(fig_dir); end

% Fixed suffixes — BS forward model
spine_suffix = '_BS';
brain_suffix = 'brain_pct10';

cmap     = make_cmap(cfg.fieldtrip_path);
cmap_hot = cmap;

fprintf('\n=== FIGURES-ONLY CONFIG ===\n')
fprintf('  Save dir:   %s\n', cfg.save_dir)
fprintf('  Steps:      %d %d %d\n', plot_step1, plot_step2, plot_step4)
if cfg.doSmooth
    fprintf('  Smoothing:  ON  (spine %d mm FWHM, brain %d mm FWHM)\n', ...
        cfg.spine_smooth_fwhm_mm, cfg.brain_smooth_fwhm_mm)
else
    fprintf('  Smoothing:  OFF\n')
end
fprintf('===========================\n\n')

%% =========================================================================
%  STEP 1 — Brain-EMG coherence figures
%% =========================================================================
if plot_step1
    fprintf('\n>>> STEP 1 figures: Brain-EMG coherence\n')

groupfile = fullfile(cfg.save_dir, ...
    sprintf('groupRes_brain_DICS_%s.mat', brain_suffix));
    assert(exist(groupfile,'file')==2, ...
        'Step 1 results not found:\n  %s', groupfile)
    load(groupfile, 'subjResults')

    % Participant 1 brain surface plot
    load(cfg.geomfile)
    mesh_brain.unit = 'mm';
    invpthr = -log10(0.05);
    sub       = subjResults(1).subjID;
    brain_int = subjResults(1).brain_int;
    mask      = subjResults(1).sig_mask > 0;
    hfig = plot_brain_surface(brain_int, mask, brain_int.coh, invpthr, cmap, ...
        mesh_brain, 'Participant 1: brain-EMG coherence', ...
        cfg.doSmooth, cfg.brain_smooth_fwhm_mm);
    if cfg.saveFigs
        savefig(hfig, fullfile(fig_dir, sprintf('step1_sub%s_brainEMG_coherence.fig', sub)))
    end

    % Group prevalence figure
    run_group_brain(subjResults, cfg.geomfile, cfg, cmap, ...
        'brain-EMG coherence', cfg.doSmooth, cfg.brain_smooth_fwhm_mm, ...
        cfg.saveFigs, fig_dir, 'step1')

    fprintf('>>> STEP 1 figures complete.\n\n')
end

%% =========================================================================
%  STEP 2 — Spinal cord-EMG coherence figures
%% =========================================================================
if plot_step2
    fprintf('\n>>> STEP 2 figures: Spinal cord-EMG coherence\n')

groupfile = fullfile(cfg.save_dir, ...
    sprintf('groupRes_spine_DICS_%s.mat', spine_suffix));

    assert(exist(groupfile,'file')==2, ...
        'Step 2 results not found:\n  %s', groupfile)
    load(groupfile, 'subjResults')
    load(cfg.geomfile)
    mesh_cut = clip_torso(mesh_torso);

    invpthr = -log10(0.05);
    nSubs   = length(subjResults);

    % subjID not stored in spine subjResults — use known subject list by index
    subs_spine = {'OP00212','OP00213','OP00215','OP00219', ...
                  'OP00220','OP00221','OP00224','OP00225','OP00226'};

    % Participant 1 spine mesh-only figure
    % (uses the group prevalence map from run_group_spine for colour accuracy;
    %  here we just plot the mesh with peak markers from stored sig_mask)
    sub      = subs_spine{1};
    mask     = subjResults(1).sig_mask;
    coh_diff    = subjResults(1).coh_diff;
    invp_smooth = subjResults(1).invp_smooth;
    pos      = subjResults(1).pos;
    thr95    = subjResults(1).thr95;
    invp_max = max(invp_smooth);
    if invp_max <= invpthr
        clim_spine = [invpthr invpthr + 0.5];
    else
        clim_spine = [invpthr invp_max];
    end

    source_p = []; source_p.pos = pos; source_p.inside = subjResults(1).inside;
    source_p.avg.coh = invp_smooth;
    cfg_interp = []; cfg_interp.parameter = 'coh';
    spine_int_p1 = ft_sourceinterpolate(cfg_interp, source_p, mesh_wm);

        if true  % participant 1 mesh-only block
            hfig_spine_mesh = figure('Color','w');
            cfg_m = []; cfg_m.figure = 'gcf'; cfg_m.method = 'surface';
            cfg_m.funparameter = 'coh'; cfg_m.funcolormap = cmap_hot;
            cfg_m.funcolorlim = clim_spine; cfg_m.projmethod = 'nearest';
            cfg_m.surffile = mesh_wm;
            ft_sourceplot(cfg_m, spine_int_p1);
            colorbar off;
            view(-250,-1); camlight; ax = gca; ax.FontSize = 14;
            title('Participant 1 — spinal cord-EMG coherence','Interpreter','none','FontSize',13);
            if cfg.saveFigs
                savefig(hfig_spine_mesh, fullfile(fig_dir, sprintf('step2_sub%s_spineEMG_meshonly.fig', sub)))
            end
        end

    % Group figures — reload sources_cent from geomfile (already loaded above)
    run_group_spine(subjResults, sources_cent, cfg.save_dir, spine_suffix, cmap_hot, ...
        mesh_wm, mesh_brain, mesh_cut, mesh_bone, mesh_lungs, mesh_heart, ...
        cfg.doSmooth, cfg.spine_smooth_fwhm_mm, cfg.saveFigs, fig_dir)

    fprintf('>>> STEP 2 figures complete.\n\n')
end


if plot_step4
    fprintf('\n>>> STEP 4 figures: Brain-spinal cord coherence boxplot\n')

    boxfile = fullfile(cfg.save_dir, ...
        'brainspine_boxplot_data_bslaw_prevalence.mat');
    assert(exist(boxfile,'file')==2, ...
        'Step 4 boxplot data not found:\n  %s', boxfile)
    load(boxfile, 'boxplot_data');

    ratio_BS     = boxplot_data.ratio_BS;
    pass_flag_BS = boxplot_data.pass_flag_BS;
    nSubs        = boxplot_data.nSubs;
    p1_idx       = boxplot_data.p1_idx;

    col_BS   = [0.8 0.2 0.4];
    col_pass = col_BS;
    col_fail = [0.65 0.65 0.65];
    col_edge = [0.15 0.15 0.15];
    col_box  = [0.93 0.93 0.93];
    col_line = [0.20 0.20 0.20];
    col_thr  = [0.10 0.10 0.10];

    q1  = quantile(ratio_BS, 0.25);
    q3  = quantile(ratio_BS, 0.75);
    med = median(ratio_BS);
    wlo = min(ratio_BS);
    whi = max(ratio_BS);

    box_hw    = 0.22;
    cap_hw    = 0.10;
    jitter_BS = linspace(-0.16, 0.16, nSubs);
    y_top     = ceil(max([ratio_BS(:); 1]) * 1.15 * 2) / 2;

    hfig_box = figure('Color','w','Position',[100 100 600 480]);
    ax = axes('Parent', hfig_box); hold(ax,'on');

    plot(ax, [0.3 1.7], [1 1], '--', 'Color', col_thr, 'LineWidth', 1.5, 'HandleVisibility','off');
    plot(ax, [1 1], [q3 whi], '-', 'Color', col_line, 'LineWidth', 1.2, 'HandleVisibility','off');
    plot(ax, [1 1], [wlo q1], '-', 'Color', col_line, 'LineWidth', 1.2, 'HandleVisibility','off');
    plot(ax, [1-cap_hw 1+cap_hw], [whi whi], '-', 'Color', col_line, 'LineWidth', 1.2, 'HandleVisibility','off');
    plot(ax, [1-cap_hw 1+cap_hw], [wlo wlo], '-', 'Color', col_line, 'LineWidth', 1.2, 'HandleVisibility','off');
    fill(ax, [1-box_hw 1+box_hw 1+box_hw 1-box_hw], [q1 q1 q3 q3], col_box, ...
        'EdgeColor', col_line, 'LineWidth', 1.2, 'HandleVisibility','off');
    plot(ax, [1-box_hw 1+box_hw], [med med], '-', 'Color', col_line, 'LineWidth', 2.2, 'HandleVisibility','off');

    for ss = 1:nSubs
        if pass_flag_BS(ss), c = col_pass; else, c = col_fail; end
        plot(ax, 1+jitter_BS(ss), ratio_BS(ss), 'o', 'MarkerFaceColor', c, ...
            'MarkerEdgeColor', col_edge, 'MarkerSize', 9, 'LineWidth', 1, 'HandleVisibility','off');
    end
    plot(ax, 1+jitter_BS(p1_idx), ratio_BS(p1_idx), 'o', 'MarkerFaceColor', 'none', ...
        'MarkerEdgeColor', 'k', 'MarkerSize', 18, 'LineWidth', 2, 'HandleVisibility','off');

    h_pass = plot(ax, nan, nan, 'o', 'MarkerFaceColor', col_pass, 'MarkerEdgeColor', col_edge, 'MarkerSize', 9, 'LineWidth', 1);
    h_fail = plot(ax, nan, nan, 'o', 'MarkerFaceColor', col_fail, 'MarkerEdgeColor', col_edge, 'MarkerSize', 9, 'LineWidth', 1);
    h_thr  = plot(ax, [NaN NaN], [NaN NaN], '--', 'Color', col_thr, 'LineWidth', 1.5);
    legend(ax, [h_pass h_fail h_thr], {'Significant','Nonsignificant','Threshold'}, ...
        'Location','best', 'FontSize', 10, 'Box','off');

    xlim(ax, [0.3 1.7]); ylim(ax, [0 y_top]);
    set(ax, 'XTick', 1, 'XTickLabel', 'Participants');
    set(ax, 'YTick', 0:0.5:y_top);
    set(ax, 'FontSize', 11, 'TickDir', 'out', 'LineWidth', 1);
    box(ax, 'off');
    ylabel(ax, 'Peak coherence / threshold', 'FontSize', 13, 'FontWeight', 'bold');
    title(ax, 'Brain-spinal cord coherence', 'FontSize', 15, 'FontWeight', 'bold', ...
        'Color', [0.15 0.15 0.15]);

    if cfg.saveFigs
        savefig(hfig_box, fullfile(fig_dir, 'brainspine_peak_vs_threshold_boxplot.fig'));
        saveas(hfig_box,  fullfile(fig_dir, 'brainspine_peak_vs_threshold_boxplot.png'));
    end
    fprintf('>>> STEP 4 figures complete.\n\n')
end

fprintf('\n=== FIGURES COMPLETE ===\n\n')


%% =========================================================================
%  GROUP ANALYSIS HELPERS  (copied verbatim from RUN_PIPELINE.m)
%% =========================================================================

function run_group_brain(subjResults, geomfile, cfg_pipeline, cmap, ...
                          label_str, doSmooth, fwhm, saveFigs, fig_dir, step_tag)

load(geomfile); mesh_brain.unit = 'mm';

nSubs     = length(subjResults);
all_masks = cat(2, subjResults(:).sig_mask);
sig_pos   = false(nSubs,1);
for ss = 1:nSubs
    if any(subjResults(ss).coh_diff > subjResults(ss).thr95)
        sig_pos(ss) = true;
    end
end

if doSmooth
    fprintf('  %g/%g subjects show sig %s (smoothed %d mm)\n', ...
        sum(sig_pos), nSubs, label_str, fwhm)
else
    fprintf('  %g/%g subjects show sig %s\n', sum(sig_pos), nSubs, label_str)
end

group_prevalence = mean(all_masks, 2);
threshold        = 0.3;
gp_masked        = group_prevalence; gp_masked(gp_masked < threshold) = 0;

group_ft     = subjResults(1).brain_int;
group_ft.pow = gp_masked;

load(cfg_pipeline.T_mat_path); T_inv = inv(T);
maxval = max(group_ft.pow); maxidx = find(group_ft.pow == maxval);
maxpos = group_ft.pos(maxidx,:);
maxpos_h = [maxpos, ones(size(maxpos,1),1)]';
x_mni = (T_inv * maxpos_h)'; x_mni = x_mni(:,1:3);
if length(maxidx)>1, disp('multiple group max locs'); end


cfg = []; cfg.parameter = 'pow'; cfg.interpmethod = 'sphere_avg'; cfg.sphereradius = 10;
group_int = ft_sourceinterpolate(cfg, group_ft, mesh_brain);

hfig_prev = figure;
cfg2 = []; cfg2.method = 'surface'; cfg2.funparameter = 'pow';
cfg2.funcolorlim = [threshold max(group_int.pow)];
cfg2.funcolormap = cmap; cfg2.projmethod = 'nearest'; cfg2.surffile = mesh_brain;
ft_sourceplot(cfg2, group_int);
view(176,-10); camlight; ax = gca; ax.FontSize = 14;
hpatch = findobj(gcf,'Type','patch'); set(hpatch,'FaceAlpha',0.9)
title(sprintf('Group prevalence — %s', label_str),'Interpreter','none')
if saveFigs, savefig(hfig_prev, fullfile(fig_dir, sprintf('%s_group_prevalence_%s.fig', step_tag, strrep(label_str,' ','_')))); end

end


function run_group_spine(subjResults, sources_cent, save_dir, out_suffix, cmap_hot, ...
                          mesh_wm, mesh_brain, mesh_cut, mesh_bone, mesh_lungs, mesh_heart, ...
                          doSmooth, fwhm, saveFigs, fig_dir)

nSubjects     = length(subjResults);
nsourcepoints = size(sources_cent.pos,1);
all_masks     = zeros(nsourcepoints, nSubjects);
for s = 1:nSubjects
    m = double(subjResults(s).coh_diff > subjResults(s).thr95);
    all_masks(:,s) = m(:);
end

sig_pos = false(nSubjects,1);
for s = 1:nSubjects
    sig_pos(s) = any(all_masks(:,s));
end

if doSmooth
    fprintf('  %d/%d subjects show sig spine-EMG coherence (smoothed %d mm)\n', ...
        sum(sig_pos), nSubjects, fwhm)
else
    fprintf('  %d/%d subjects show sig spine-EMG coherence\n', sum(sig_pos), nSubjects)
end

threshold      = 0.2;
prevalence_loc = mean(all_masks, 2);
group_ft = []; group_ft.pos = sources_cent.pos;
group_ft.inside = sources_cent.inside;
group_ft.pow = prevalence_loc; group_ft.pow(group_ft.pow < threshold) = 0;

cfg = []; cfg.parameter = 'pow'; cfg.interpmethod = 'nearest';
group_int = ft_sourceinterpolate(cfg, group_ft, mesh_wm);

hfig_spineprev = figure;
cfg2 = []; cfg2.method = 'surface'; cfg2.funparameter = 'pow';
cfg2.maskparameter = 'mask'; cfg2.funcolorlim = [threshold max(group_int.pow)];
cfg2.funcolormap = cmap_hot; cfg2.projmethod = 'nearest'; cfg2.surffile = mesh_wm;
cfg2.opacitylim = [threshold max(group_int.pow)]; cfg2.opacitymap = 'rampup';
ft_sourceplot(cfg2, group_int);
view(-250,-1); camlight; ax = gca; ax.FontSize = 14; hold on
ft_plot_mesh(mesh_brain,'facecolor',[0.8 0.3 0.3],'facealpha',0.07,'edgecolor','none');
ft_plot_mesh(mesh_cut,  'facecolor',[0.3 0.3 0.9],'facealpha',0.1, 'edgecolor','none');
ft_plot_mesh(mesh_bone, 'facecolor',[0.9 0.85 0.7],'facealpha',0.3, 'edgecolor','none');
ft_plot_mesh(mesh_lungs,'facecolor',[0.8 0.3 0.3],'facealpha',0.1, 'edgecolor','none');
ft_plot_mesh(mesh_heart,'facecolor',[0.8 0.3 0.3],'facealpha',0.1, 'edgecolor','none');
title('Group prevalence — spinal cord-EMG','Interpreter','none')
if saveFigs, savefig(hfig_spineprev, fullfile(fig_dir, 'step2_group_spineEMG_prevalence.fig')); end

% Group prevalence — mesh only, no colorbar
hfig_grp_mesh = figure('Color','w','Position',[100 100 400 650]);
cfg3 = []; cfg3.method = 'surface'; cfg3.funparameter = 'pow';
cfg3.funcolorlim = [threshold max(group_int.pow)];
cfg3.funcolormap = cmap_hot; cfg3.projmethod = 'nearest'; cfg3.surffile = mesh_wm;
cfg3.opacitylim = [threshold max(group_int.pow)]; cfg3.opacitymap = 'rampup';
ft_sourceplot(cfg3, group_int);
colorbar off;
view(-250,-1); camlight; ax = gca; ax.FontSize = 14;
title('Group prevalence — spinal cord-EMG','Interpreter','none','FontSize',13)
if saveFigs, savefig(hfig_grp_mesh, fullfile(fig_dir, 'step2_group_spineEMG_prevalence_meshonly.fig')); end


end


%% =========================================================================
%  SHARED UTILITY FUNCTIONS  (copied verbatim from RUN_PIPELINE.m)
%% =========================================================================

function hfig = plot_brain_surface(brain_int, mask, invp, invpthr, cmap, mesh_brain, ttl, doSmooth, fwhm)
    hfig = figure;
    cfg = []; cfg.figure = 'gcf'; cfg.method = 'surface';
    cfg.funparameter = 'coh'; cfg.funcolormap = cmap;
    if any(mask)
        clim_upper = max(invp(mask));
        if clim_upper <= invpthr
            clim_upper = invpthr + 0.1;
        end
        cfg.funcolorlim = [invpthr clim_upper];
    else
        cfg.funcolorlim = [invpthr invpthr+1];
    end
    cfg.projmethod = 'nearest'; cfg.surffile = mesh_brain;
    ft_sourceplot(cfg, brain_int);
    view(176,-10); camlight; ax = gca; ax.FontSize = 14;
    hpatch = findobj(gcf,'Type','patch'); set(hpatch,'FaceAlpha',0.9)
    title(ttl,'Interpreter','none')
end

function mesh_cut = clip_torso(mesh_torso)
    y = mesh_torso.vertices(:,2);
    keep_vert = y > -200;
    new_idx = zeros(size(keep_vert)); new_idx(keep_vert) = 1:sum(keep_vert);
    faces_keep = all(keep_vert(mesh_torso.faces), 2);
    mesh_cut.vertices = mesh_torso.vertices(keep_vert,:);
    mesh_cut.faces    = new_idx(mesh_torso.faces(faces_keep,:));
    mesh_cut.unit     = mesh_torso.unit;
end

function cmap = make_cmap(fieldtrip_path)
    ncol = 256;
    addpath(fullfile(fieldtrip_path,'external','matplotlib'))
    cmap = [[0.92 0.92 0.92]; flipud(magma(ncol-1))];
end