%% step3_VE_BS_prevalenceROI.m
% Step 3: Build spinal virtual electrode using the Biot-Savart (BS)
% leadfield and the data-driven prevalence-cluster ROI from step 2
% for all subjects.
%
% Inputs (must exist before running):
%   geomfile           - geometries_experimental.mat: sources_cent, mesh_torso, mesh_wm
%   lf_path            - leadfield_experimental_bslaw_experimental.mat: leadfield_bs
%   cluster_path       - cluster_spineEMG_pos_BS.mat: ROIpos (from step 2 group analysis)
%   data_root          - per-subject SPM .mat files (see datafile pattern below)
%
% Outputs (written to save_dir / fig_dir):
%   VE_spine_prevalence_sub<ID>_forspectra_BS.mat    - per-subject VE, roi_center, R, ROIpos
%   step3_prevalenceROI_VE_BS.fig / .png            - ROI location on mesh

clear all; close all; clc;

%% =========================================================================
%  USER CONFIG
%% =========================================================================
% Machine-specific paths live in source/brainspine_config.m — edit that
% file to match your local installation.
repo_root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(repo_root, 'source'));
paths = brainspine_config();

fieldtrip_path = paths.fieldtrip_path;
spm_path       = paths.spm_path;
bsc_path       = fullfile(repo_root, 'source');
data_root      = paths.data_root;
save_dir       = paths.save_dir;
geomfile       = paths.geomfile;
lf_path        = paths.lf_path;
cluster_path   = fullfile(save_dir, 'cluster_spineEMG_pos_BS.mat');

subs_spine = {'OP00212','OP00213','OP00215','OP00219', ...
              'OP00220','OP00221','OP00224','OP00225','OP00226'};

out_suffix = '_BS';

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
%  LOAD GEOMETRY AND PREVALENCE-CLUSTER ROI
%% =========================================================================
fprintf('Loading geometry...\n');
geom_exp     = load(geomfile);
sources_cent = geom_exp.sources_cent;
mesh_torso   = geom_exp.mesh_torso;
mesh_wm      = geom_exp.mesh_wm;
mesh_wm.unit = 'mm';

fprintf('Loading prevalence-cluster ROI from step 2...\n');
cluster_data = load(cluster_path, 'ROIpos');
ROIpos       = cluster_data.ROIpos;

roi_center = mean(ROIpos, 1);
R          = max(sqrt(sum((ROIpos - roi_center).^2, 2)));

fprintf('  Prevalence-cluster ROI: %d points\n', size(ROIpos,1));
fprintf('  ROI centre: [%.1f %.1f %.1f] mm\n', roi_center);
fprintf('  ROI radius: %.1f mm\n', R);

% Plot ROI on mesh
hfig_roi = figure('Color','w');
ft_plot_mesh(mesh_wm,'facecolor',[0.7 0.7 0.7],'facealpha',0.3,'edgecolor','none');
hold on;
plot3(ROIpos(:,1), ROIpos(:,2), ROIpos(:,3), 'o', 'MarkerSize', 10, ...
    'MarkerEdgeColor',[0.9 0.3 0], 'MarkerFaceColor',[1 0.4 0.1], 'LineWidth', 2);
view(90,18); material dull;
title('Spinal ROI from prevalence cluster (BS)','Interpreter','none');
savefig(hfig_roi, fullfile(fig_dir, ['step3_prevalenceROI_VE' out_suffix '.fig']));
saveas(hfig_roi,  fullfile(fig_dir, ['step3_prevalenceROI_VE' out_suffix '.png']));
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

    %% Build leadfield — label-based matching (BS)
    fprintf('  Building leadfield...\n');
    lf_data = load(lf_path);
    lf_raw  = lf_data.leadfield_bs;

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

    %% Extract virtual electrode from prevalence-cluster ROI
    cfg_ve = [];
    cfg_ve.pos          = roi_center;
    cfg_ve.radius       = R;
    cfg_ve.method       = 'svd';
    cfg_ve.numcomponent = 1;
    VE = ft_virtualchannel_sphere(cfg_ve, ftdat, source_idx);

    savename = sprintf('VE_spine_prevalence_sub%s_forspectra%s', sub, out_suffix);
    save(fullfile(save_dir, savename), 'VE', 'roi_center', 'R', 'ROIpos');
    fprintf('  VE saved: %s\n', savename);
end

fprintf('\n=== STEP 3 DONE ===\n');
