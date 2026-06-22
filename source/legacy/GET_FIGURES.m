%% RUN_FIGURES_ONLY.m
% Regenerates all pipeline figures from pre-computed results.
%
% Requires the following saved files in save_dir:
%   groupRes_brain_DICS_bemv2_brainEMG_brainSmooth_8mm.mat
%   groupRes_spine_DICS_bemv2_permSmooth_20mm.mat
%   cluster_spineEMG_pos_bemv2_permSmooth_20mm.mat
%   groupRes_brain_DICS_spineVC_bemv2_functionalVE_spineSmooth_20mm_brainSmooth_8mm.mat
%   M1_ROI_bemv2_functionalVE_spineSmooth_20mm_brainSmooth_8mm.mat
%   geometries_cervical_realistic.mat  (geomfile)
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
cfg.geomfile        = 'C:\Users\mspedden\Documents\Leadfields meshes\geometries_cervical_realistic.mat';
cfg.T_mat_path      = 'C:\Users\mspedden\Documents\Leadfields meshes\T.mat';
cfg.save_dir        = 'C:\Users\mspedden\Documents\forGeneratingFigures';

% --- Smoothing config (must match what was used to generate the results) --
cfg.doSmooth             = 1;
cfg.spine_smooth_fwhm_mm = 20;
cfg.brain_smooth_fwhm_mm = 8;


% --- Steps to plot --------------------------------------------------------
plot_step1 = 1;   % Brain-EMG coherence figures
plot_step2 = 1;   % Spine-EMG coherence figures
plot_step3 = 1;   % Spinal VE ROI figure
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

% Filename suffixes — derived from smoothing config
if cfg.doSmooth
    spine_suffix  = sprintf('_permSmooth_%dmm', cfg.spine_smooth_fwhm_mm);
    brain_suffix  = sprintf('brainEMG_brainSmooth_%dmm', cfg.brain_smooth_fwhm_mm);
    ve_suffix     = sprintf('functionalVE_spineSmooth_%dmm_brainSmooth_%dmm', ...
                        cfg.spine_smooth_fwhm_mm, cfg.brain_smooth_fwhm_mm);
else
    spine_suffix  = '';
    brain_suffix  = 'brainEMG';
    ve_suffix     = 'functionalVE';
end

cmap     = make_cmap(cfg.fieldtrip_path);
cmap_hot = cmap;

fprintf('\n=== FIGURES-ONLY CONFIG ===\n')
fprintf('  Save dir:   %s\n', cfg.save_dir)
fprintf('  Steps:      %d %d %d %d\n', plot_step1, plot_step2, plot_step3, plot_step4)
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
        sprintf('groupRes_brain_DICS_bemv2_%s.mat', brain_suffix));
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
%  STEP 2 — Spine-EMG coherence figures
%% =========================================================================
if plot_step2
    fprintf('\n>>> STEP 2 figures: Spine-EMG coherence\n')

    groupfile = fullfile(cfg.save_dir, ...
        sprintf('groupRes_spine_DICS_bemv2%s.mat', spine_suffix));
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
    clim_spine = [2.3 2.47];

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
            title('Participant 1 — spine-EMG coherence','Interpreter','none','FontSize',13);
            peak_tol = 0.01;
            sig_invp = invp_smooth; sig_invp(~mask) = -inf;
            peaks_p1 = pos(sig_invp >= max(sig_invp(isfinite(sig_invp)))*(1-peak_tol), :);
            hold on;
            scatter3(peaks_p1(:,1), peaks_p1(:,2), peaks_p1(:,3)+10, ...
                200, '.', 'filled', 'MarkerFaceColor',[1 1 0], 'MarkerEdgeColor','k','LineWidth',1.5);
            scatter_obj = findobj(gca,'Type','Scatter'); uistack(scatter_obj,'top');
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

%% =========================================================================
%  STEP 3 — Spinal VE ROI figure
%% =========================================================================
if plot_step3
    fprintf('\n>>> STEP 3 figures: Spinal VE ROI\n')

    roi_file = fullfile(cfg.save_dir, ...
        ['cluster_spineEMG_pos_bemv2' spine_suffix '.mat']);
    assert(exist(roi_file,'file')==2, ...
        'Step 3 ROI file not found:\n  %s', roi_file)
    load(roi_file)   % loads ROIpos
    load(cfg.geomfile)

    hfig_roi = figure('Color','w');
    ft_plot_mesh(mesh_wm,'facecolor',[0.7 0.7 0.7],'facealpha',0.3,'edgecolor','none');
    hold on
    plot3(ROIpos(:,1),ROIpos(:,2),ROIpos(:,3),'o','MarkerSize',10, ...
        'MarkerEdgeColor',[0.9 0.3 0],'MarkerFaceColor',[1 0.4 0.1],'LineWidth',2);
    view(90,18); material dull;
    if cfg.doSmooth
        title(sprintf('Spinal ROI (smoothed %d mm)', cfg.spine_smooth_fwhm_mm),'Interpreter','none')
    else
        title('Spinal ROI','Interpreter','none')
    end
    if cfg.saveFigs, savefig(hfig_roi, fullfile(fig_dir, 'step3_spinal_VE_ROI.fig')); end

    fprintf('>>> STEP 3 figures complete.\n\n')
end

%% =========================================================================
%  STEP 4 — Brain-SpineVE coherence + M1 overlap figures
%% =========================================================================
if plot_step4
    fprintf('\n>>> STEP 4 figures: Brain-SpineVE coherence\n')

    groupfile = fullfile(cfg.save_dir, ...
        sprintf('groupRes_brain_DICS_spineVC_bemv2_%s.mat', ve_suffix));
    assert(exist(groupfile,'file')==2, ...
        'Step 4 results not found:\n  %s', groupfile)
    load(groupfile, 'subjResults')

    load(cfg.geomfile)
    mesh_brain.unit = 'mm';
    invpthr = -log10(0.05);

    % Participant 1 brain surface plot
    sub       = subjResults(1).subjID;
    brain_int = subjResults(1).brain_int;
    mask      = subjResults(1).sig_mask > 0;
    hfig = plot_brain_surface(brain_int, mask, brain_int.coh, invpthr, cmap, ...
        mesh_brain, 'Participant 1: brain-spineVE coherence', ...
        cfg.doSmooth, cfg.brain_smooth_fwhm_mm);
    if cfg.saveFigs
        savefig(hfig, fullfile(fig_dir, sprintf('step4_sub%s_brainSpineVE_coherence.fig', sub)))
    end

    % Group prevalence figure
    run_group_brain(subjResults, cfg.geomfile, cfg, cmap, ...
        'brain-spineVE coherence', cfg.doSmooth, cfg.brain_smooth_fwhm_mm, ...
        cfg.saveFigs, fig_dir, 'step4')


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

group_ft = []; group_ft.pos = subjResults(1).pos;
group_ft.inside = subjResults(1).inside; group_ft.pow = gp_masked;

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
if doSmooth
    title(sprintf('Group prevalence — %s (smoothed %d mm)', label_str, fwhm),'Interpreter','none')
else
    title(sprintf('Group prevalence — %s', label_str),'Interpreter','none')
end
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
if doSmooth
    title(sprintf('Group prevalence — spine-EMG (smoothed %d mm)', fwhm),'Interpreter','none')
else
    title('Group prevalence — spine-EMG','Interpreter','none')
end
if saveFigs, savefig(hfig_spineprev, fullfile(fig_dir, 'step2_group_spineEMG_prevalence.fig')); end

% Group prevalence — mesh only, no colorbar, star at peak(s)
hfig_grp_mesh = figure('Color','w','Position',[100 100 400 650]);
cfg3 = []; cfg3.method = 'surface'; cfg3.funparameter = 'pow';
cfg3.funcolorlim = [threshold max(group_int.pow)];
cfg3.funcolormap = cmap_hot; cfg3.projmethod = 'nearest'; cfg3.surffile = mesh_wm;
cfg3.opacitylim = [threshold max(group_int.pow)]; cfg3.opacitymap = 'rampup';
ft_sourceplot(cfg3, group_int);
colorbar off;
view(-250,-1); camlight; ax = gca; ax.FontSize = 14;
if doSmooth
    title(sprintf('Group prevalence — spine-EMG (smoothed %d mm)', fwhm),'Interpreter','none','FontSize',13)
else
    title('Group prevalence — spine-EMG','Interpreter','none','FontSize',13)
end
peak_tol  = 0.01;
prev_vals = group_ft.pow; prev_vals(prev_vals < threshold) = -inf;
peaks_grp = sources_cent.pos(prev_vals >= max(prev_vals(isfinite(prev_vals)))*(1-peak_tol), :);
hold on;
scatter3(peaks_grp(:,1), peaks_grp(:,2), peaks_grp(:,3)+10, ...
    200, 'p', 'filled', 'MarkerFaceColor',[1 1 0], 'MarkerEdgeColor','k','LineWidth',1.5);
scatter_obj = findobj(gca,'Type','Scatter'); uistack(scatter_obj,'top');
if saveFigs, savefig(hfig_grp_mesh, fullfile(fig_dir, 'step2_group_spineEMG_prevalence_meshonly.fig')); end

% Subject line plot
subj_cmap = [27,158,119; 217,95,2; 117,112,179; 231,41,138; 102,166,30; ...
             230,171,2;  166,118,29; 102,102,102; 55,126,184] / 255;
x = sources_cent.pos(:,2);
hfig_lines = figure; hold on;
for s = 1:nSubjects
    cdiff = subjResults(s).coh_diff;
    thr   = subjResults(s).thr95;
    sig   = cdiff > thr;
    c     = subj_cmap(s,:);
    for i = 1:length(x)-1
        if sig(i)&&sig(i+1)
            plot(x(i:i+1),cdiff(i:i+1),'-','Color',c,'LineWidth',1.5,'HandleVisibility','off')
        else
            plot(x(i:i+1),cdiff(i:i+1),'-','Color',[0.7 0.7 0.7],'HandleVisibility','off')
        end
    end
    plot(x(sig),cdiff(sig),'.','Color',c,'MarkerSize',12,'HandleVisibility','off')
    if sig_pos(s)
        h(s) = plot(nan,nan,'-','Color',c,'LineWidth',1.5);
    else
        h(s) = plot(nan,nan,'-','Color',[0.7 0.7 0.7],'LineWidth',1.5);
    end
end
yline(0,':k','HandleVisibility','off')
xlabel('Cranial-caudal position (mm)'); ylabel('Coherence difference');
if doSmooth
    title(sprintf('Spine-EMG coherence differences (smoothed %d mm)', fwhm))
else
    title('Spine-EMG coherence differences')
end
legend(h, arrayfun(@(s) sprintf('Participant %d',s), 1:nSubjects,'UniformOutput',false), ...
    'Location','bestoutside');
set(gca,'FontSize',13); grid on;
if saveFigs, savefig(hfig_lines, fullfile(fig_dir, 'step2_group_spineEMG_subject_lines.fig')); end


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
    if doSmooth
        title(sprintf('%s (smoothed %d mm)', ttl, fwhm),'Interpreter','none')
    else
        title(ttl,'Interpreter','none')
    end
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