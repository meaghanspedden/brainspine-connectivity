%% step_3_VE_true_bem_anatROI.m
% Step 3: Build spinal virtual electrode from anatomically defined ROI
% (sources 14:27, C6-T1 spinal segments) for all subjects.

clear all; close all; clc;

%% =========================================================================
%  USER CONFIG
%% =========================================================================
fieldtrip_path = 'C:\Users\mspedden\Documents\fieldtrip';
spm_path       = 'C:\Users\mspedden\Documents\spm';
bsc_path       = 'C:\Users\mspedden\Documents\brainspineconnectivity\source';
data_root      = 'C:\spinecoh_data';
save_dir       = 'C:\Users\mspedden\Documents\brainspine_savetest';
geomfile       = 'C:\Leadfields meshes\geometries_experimental.mat';
lf_path        = 'C:\Leadfields meshes\leadfield_experimental_bem_experimental.mat';

subs_spine = {'OP00212','OP00213','OP00215','OP00219', ...
              'OP00220','OP00221','OP00224','OP00225','OP00226'};

out_suffix = '_true_bem';
roi_idx    = 14:27;   % C6-T1 spinal segments (anatomically defined)

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
geom_exp     = load(geomfile);
sources_cent = geom_exp.sources_cent;
mesh_torso   = geom_exp.mesh_torso;
mesh_wm      = geom_exp.mesh_wm;
mesh_wm.unit = 'mm';

% ROI positions and centre
ROIpos     = sources_cent.pos(roi_idx, :);
roi_center = mean(ROIpos, 1);
R          = max(sqrt(sum((ROIpos - roi_center).^2, 2)));

fprintf('  Anatomical ROI: sources %d-%d (%d points)\n', ...
    roi_idx(1), roi_idx(end), numel(roi_idx));
fprintf('  ROI centre: [%.1f %.1f %.1f] mm\n', roi_center);
fprintf('  ROI radius: %.1f mm\n', R);

% Plot ROI on mesh
hfig_roi = figure('Color','w');
ft_plot_mesh(mesh_wm,'facecolor',[0.7 0.7 0.7],'facealpha',0.3,'edgecolor','none');
hold on;
plot3(ROIpos(:,1), ROIpos(:,2), ROIpos(:,3), 'o', 'MarkerSize', 10, ...
    'MarkerEdgeColor',[0.9 0.3 0], 'MarkerFaceColor',[1 0.4 0.1], 'LineWidth', 2);
view(90,18); material dull;
title('Anatomical spinal ROI C6-T1 (true BEM)','Interpreter','none');
savefig(hfig_roi, fullfile(fig_dir, ['step3_anatROI_VE' out_suffix '.fig']));
saveas(hfig_roi,  fullfile(fig_dir, ['step3_anatROI_VE' out_suffix '.png']));
close(hfig_roi);

%% =========================================================================
%  SUBJECT LOOP
%% =========================================================================
for ss = 1:length(subs_spine)
    sub = subs_spine{ss};
    fprintf('\n========================================\n');
    fprintf('  Subject %s (%d/%d)\n', sub, ss, length(subs_spine));
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

    % Rectify EMG
    cfg = []; cfg.rectify = 'yes'; cfg.channel = 'EXG1';
    ftdatr = ft_preprocessing(cfg, ftdat);
    for k = 1:length(ftdat.trial)
        ftdat.trial{k}(end,:) = ftdatr.trial{k};
    end

    %% Build leadfield — label-based matching
    fprintf('  Building leadfield...\n');
    lf_data = load(lf_path);
    lf_raw  = lf_data.leadfield_cord;

    data_meg_labels        = ftdat.label(~strcmp(ftdat.label,'EXG1'));
    [common_labels,idx_lf] = intersect(lf_raw.label, data_meg_labels, 'stable');
    fprintf('  Data MEG: %d  |  LF: %d  |  Matched: %d\n', ...
        numel(data_meg_labels), numel(lf_raw.label), numel(common_labels));

    Lf        = lf_raw;
    Lf.label  = common_labels;
    Lf.pos    = sources_cent.pos;
    Lf.inside = ones(size(sources_cent.pos,1), 1);
    for i = 1:numel(lf_raw.leadfield)
        if ~isempty(lf_raw.leadfield{i})
            Lf.leadfield{i} = lf_raw.leadfield{i}(idx_lf, :);
        end
    end

    %% Volume conductor
    cfg = []; cfg.method = 'infinite'; cfg.siunits = 1;
    cfg.grad = grad_mm; cfg.conductivity = 1;
    dummyvol = ft_prepare_headmodel(cfg, mesh_torso);

    %% LCMV beamformer
    cfg = []; cfg.covariance = 'yes';
    cfg.covariancewindow = 'all';
    cfg.channel = Lf.label;
    tlock = ft_timelockanalysis(cfg, ftdat);

    cfg_lcmv = [];
    cfg_lcmv.method                   = 'lcmv';
    cfg_lcmv.headmodel                = dummyvol;
    cfg_lcmv.sourcemodel.pos          = Lf.pos;
    cfg_lcmv.sourcemodel.unit         = 'mm';
    cfg_lcmv.sourcemodel.inside       = logical(Lf.inside);
    cfg_lcmv.sourcemodel.leadfield    = Lf.leadfield;
    cfg_lcmv.sourcemodel.label        = Lf.label;
    cfg_lcmv.lcmv.keepfilter          = 'yes';
    source_idx = ft_sourceanalysis(cfg_lcmv, tlock);

    %% Extract virtual electrode from anatomical ROI
    cfg_ve = [];
    cfg_ve.pos          = roi_center;
    cfg_ve.radius       = R;
    cfg_ve.method       = 'svd';
    cfg_ve.numcomponent = 1;
    VE = ft_virtualchannel(cfg_ve, ftdat, source_idx);

    savename = sprintf('VE_spine_anatROI_sub%s_forspectra_bemv2%s', sub, out_suffix);
    save(fullfile(save_dir, savename), 'VE', 'roi_center', 'R', 'ROIpos');
    fprintf('  VE saved: %s\n', savename);
end

fprintf('\n=== STEP 3 DONE ===\n');