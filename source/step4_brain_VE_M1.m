%% step4_brain_VE_M1.m
% Step 4: Build brain virtual electrode from M1 ROI.
% M1 ROI = intersection of brain-EMG significant region and M1 sphere.
%
% Depends on step1's output (groupRes_brain_DICS_brain_pct10.mat) for the
% group prevalence map used to define the M1 ROI.

clear all; close all; clc;

%% =========================================================================
%  USER CONFIG
%% =========================================================================
fieldtrip_path = 'C:\Users\mspedden\Documents\fieldtrip';
spm_path       = 'C:\Users\mspedden\Documents\spm';
bsc_path       = 'C:\Users\mspedden\Documents\brainspineconnectivity\source';
data_root      = 'C:\spinecoh_data';
save_dir       = 'C:\Users\mspedden\Documents\brainspine_savetest';
geomfile       = 'C:\Leadfields meshes\geometries_experimental_withbrain.mat';
T_mat_path     = 'C:\Leadfields meshes\T.mat';

subs_brain = {'OP00212','OP00213','OP00215','OP00219', ...
              'OP00225','OP00221','OP00224'};

out_suffix     = '_brain_pct10';
M1_mni_centre  = [-38 -26 56];   % left M1 hand knob (MNI mm)
M1_radius_mm   = 20;

%% =========================================================================
%  SETUP
%% =========================================================================
addpath(bsc_path); addpath(spm_path);
spm('defaults','EEG');
addpath(fieldtrip_path);
ft_defaults;

fig_dir = fullfile(save_dir, 'figures');
if ~exist(fig_dir,'dir'), mkdir(fig_dir); end

%% =========================================================================
%  LOAD GEOMETRY
%% =========================================================================
fprintf('Loading geometry...\n');
load(geomfile);   % sources_brain, mesh_brain etc.
mesh_brain.unit = 'mm';
load(T_mat_path); T_inv = inv(T);

%% =========================================================================
%  COMPUTE M1 ROI FROM BRAIN-EMG SIGNIFICANT REGION
%% =========================================================================
fprintf('Computing M1 ROI...\n');

% Load group brain-EMG results
group_file = fullfile(save_dir, ['groupRes_brain_DICS' out_suffix '.mat']);
loaded     = load(group_file);
subjResults_brain = loaded.subjResults;

% Remove empty
valid = ~cellfun(@isempty, {subjResults_brain.subjID});
subjResults_brain = subjResults_brain(valid);
nSubs = numel(subjResults_brain);

% Group prevalence of significant sources
all_masks        = cat(2, subjResults_brain(:).sig_mask);
group_prevalence = mean(all_masks, 2);
threshold        = 0.3;
active_idx       = find(group_prevalence >= threshold);
voxel_pos        = subjResults_brain(1).brain_int.pos;
active_pos       = voxel_pos(active_idx, :);
active_pos       = active_pos(~any(isnan(active_pos),2), :);

% M1 sphere in native space
[xg,yg,zg]  = ndgrid(-M1_radius_mm:M1_radius_mm, ...
                      -M1_radius_mm:M1_radius_mm, ...
                      -M1_radius_mm:M1_radius_mm);
sphere_mask  = (xg.^2 + yg.^2 + zg.^2) <= M1_radius_mm^2;
roi_coords   = [xg(sphere_mask), yg(sphere_mask), zg(sphere_mask)] + M1_mni_centre;
roi_hom      = [roi_coords, ones(size(roi_coords,1),1)];
roi_native   = (T * roi_hom')'; roi_native = roi_native(:,1:3);

% Intersection — active voxels within M1 sphere
D_mat          = pdist2(roi_native, active_pos);
min_dist       = min(D_mat, [], 2);
roi_overlap    = min_dist <= 5;
intersection_pos = roi_native(roi_overlap, :);

fprintf('  M1 sphere: %d voxels\n', size(roi_native,1));
fprintf('  Active (prevalence>=%.1f): %d voxels\n', threshold, size(active_pos,1));
fprintf('  Intersection: %d voxels\n', size(intersection_pos,1));

if isempty(intersection_pos)
    warning('No intersection found — using full M1 sphere instead');
    intersection_pos = roi_native;
end

% Save M1 ROI
save(fullfile(save_dir, ['M1_ROI' out_suffix '.mat']), ...
    'intersection_pos', 'roi_native', 'M1_mni_centre', 'M1_radius_mm');
fprintf('  M1 ROI saved.\n');

% Plot M1 ROI on brain mesh
hfig_m1 = figure('Color','w');
ft_plot_mesh(mesh_brain,'facecolor',[0.85 0.85 0.85],'facealpha',0.3,'edgecolor','none');
hold on;
plot3(roi_native(:,1), roi_native(:,2), roi_native(:,3), 'o', ...
    'MarkerSize',4,'MarkerFaceColor',[0.6 0.6 0.6],'MarkerEdgeColor','none');
plot3(active_pos(:,1), active_pos(:,2), active_pos(:,3), 'o', ...
    'MarkerSize',5,'MarkerFaceColor',[0.2 0.4 0.8],'MarkerEdgeColor','none');
if ~isempty(intersection_pos)
    plot3(intersection_pos(:,1), intersection_pos(:,2), intersection_pos(:,3), 'o', ...
        'MarkerSize',7,'MarkerFaceColor',[0.1 0.8 0.2],'MarkerEdgeColor','k','LineWidth',0.5);
end
axis equal; camlight; lighting gouraud; view(176,-10);
legend({'M1 sphere','Sig brain-EMG','Intersection'},'Location','best','FontSize',11,'Box','off');
title(sprintf('M1 ROI — brain-EMG intersection  |  %d voxels', ...
    size(intersection_pos,1)),'Interpreter','none','FontSize',12);
savefig(hfig_m1, fullfile(fig_dir, ['step4_M1_ROI' out_suffix '.fig']));
saveas(hfig_m1,  fullfile(fig_dir, ['step4_M1_ROI' out_suffix '.png']));
close(hfig_m1);

%% =========================================================================
%  HEAD MODEL (same for all subjects — built once)
%% =========================================================================
cfg = []; cfg.method = 'singleshell';
vol = ft_prepare_headmodel(cfg, mesh_brain);

%% =========================================================================
%  SUBJECT LOOP — build brain VE from M1 ROI
%% =========================================================================
for ss = 1:length(subs_brain)
    sub = subs_brain{ss};
    fprintf('\n========================================\n');
    fprintf('  Subject %s (%d/%d)\n', sub, ss, length(subs_brain));
    fprintf('========================================\n');

    %% Load data
    run = '001'; if strcmp(sub,'OP00224'), run = '002'; end
    datafile = fullfile(data_root, ['sub-' sub], 'ses-001', 'meg', ...
        sprintf('pmergedoe1000mspddfflo45hi45hfcstatic_%s_array1.mat', run));

    D       = spm_eeg_load(datafile);
    grad_mm = D.sensors('MEG');
    ftdat   = spm2fieldtrip(D);

    badchans = D.chanlabels(D.badchannels);
    cfg = []; cfg.channel = setdiff(ftdat.label, badchans);
    ftdat = ft_selectdata(cfg, ftdat);

    % Brain sensors only
    brainidx  = find(grad_mm.chanpos(:,2) > 200);
    braingrad = subset_grad(grad_mm, brainidx);

    %% Leadfield (subject-specific grad, shared head model)
    cfg = []; cfg.sourcemodel = sources_brain; cfg.headmodel = vol;
    cfg.grad = braingrad; cfg.reducerank = 'no';
    LF = ft_prepare_leadfield(cfg);

    %% LCMV beamformer
    cfg = []; cfg.covariance = 'yes';
    tlock = ft_timelockanalysis(cfg, ftdat);

    cfg_lcmv = []; cfg_lcmv.method = 'lcmv';
    cfg_lcmv.headmodel             = vol;
    cfg_lcmv.sourcemodel.leadfield = LF;
    cfg_lcmv.grid                  = sources_brain;
    cfg_lcmv.unit                  = LF.unit;
    cfg_lcmv.lcmv.keepfilter       = 'yes';
    source_idx = ft_sourceanalysis(cfg_lcmv, tlock);

    %% Map M1 ROI to source grid
    idx_roi    = dsearchn(source_idx.pos, intersection_pos);
    idx_roi    = unique(idx_roi);
    idx_roi    = idx_roi(source_idx.inside(idx_roi));

    roi_center = mean(source_idx.pos(idx_roi,:), 1);
    d          = sqrt(sum((intersection_pos - roi_center).^2, 2));
    R          = max(d);

    fprintf('  M1 ROI: %d source points  centre=[%.1f %.1f %.1f]  R=%.1f mm\n', ...
        numel(idx_roi), roi_center, R);

    %% Extract brain VE
    cfg_ve = [];
    cfg_ve.pos          = roi_center;
    cfg_ve.radius       = R;
    cfg_ve.method       = 'svd';
    cfg_ve.numcomponent = 1;
    VE_brain = ft_virtualchannel_sphere(cfg_ve, ftdat, source_idx);

    savename = sprintf('sub%s_VE_brain_M1%s', sub, out_suffix);
    save(fullfile(save_dir, savename), 'VE_brain', 'roi_center', 'R', 'idx_roi');
    fprintf('  Brain VE saved: %s\n', savename);
end

fprintf('\n=== STEP 4 DONE ===\n');

%% =========================================================================
%  LOCAL FUNCTIONS
%% =========================================================================
function braingrad = subset_grad(grad_mm, idx)
    braingrad          = grad_mm;
    braingrad.chanori  = grad_mm.chanori(idx,:);
    braingrad.chanpos  = grad_mm.chanpos(idx,:);
    braingrad.chantype = grad_mm.chantype(1:length(idx));
    braingrad.chanunit = grad_mm.chanunit(1:length(idx));
    braingrad.coilori  = grad_mm.coilori(idx,:);
    braingrad.coilpos  = grad_mm.coilpos(idx,:);
    braingrad.label    = grad_mm.label(idx);
    braingrad.tra      = grad_mm.tra(idx,idx);
end
