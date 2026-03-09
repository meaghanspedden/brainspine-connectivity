% virtual electrode spinal cord (from prevalence analysis) — BEM v2
% Option to use ROI from unsmoothed or spatially smoothed spine analysis

clear all
close all
clc

addpath('C:\Users\mspedden\Documents\brainspineconnectivity\source')
addpath('C:\Users\mspedden\Documents\spm')
spm('defaults','EEG')

addpath('C:\Users\mspedden\Documents\fieldtrip')
ft_defaults

save_dir='C:\Users\mspedden\Documents\brainspine_save_bemv2';
if ~exist(save_dir,'dir')
    mkdir(save_dir)
end

subs = {'OP00212','OP00213','OP00215', 'OP00219', ...
    'OP00225', 'OP00221', 'OP00224'};

generic_dir = 'C:\Users\mspedden\Documents\new_leadfields_and_geom';
geomfile = fullfile(generic_dir, 'geometries_cervical_realistic.mat');

LFop    = 'spine';
HFC     = 1;
rectify = 1;
fband   = [10 35];

%% option: use ROI from smoothed or unsmoothed spine analysis
useSmoothROI = 1;          % 0 = unsmoothed ROI, 1 = smoothed ROI
smooth_fwhm_mm = 40;       % must match the smoothing used in the analysis script

if useSmoothROI
    roi_file = fullfile(save_dir, sprintf('cluster_spineEMG_pos_bemv2_permSmooth_%dmm.mat', smooth_fwhm_mm));
    ve_suffix = sprintf('_permSmooth_%dmm', smooth_fwhm_mm);
else
    roi_file = fullfile(save_dir, 'cluster_spineEMG_pos_bemv2.mat');
    ve_suffix = '';
end

% load cluster position from chosen v2 analysis
load(roi_file)  % loads ROIpos

for ss = 1:length(subs)

    sub = subs{ss};

    if strcmp(sub, 'OP00224')
        datwithEMGmerged = fullfile('C:\Users\mspedden\Documents', ...
            ['sub-' sub], 'ses-001', 'meg', ...
            'pmergedoe1000mspddfflo45hi45hfcstatic_002_array1.mat');
    else
        datwithEMGmerged = fullfile('C:\Users\mspedden\Documents', ...
            ['sub-' sub], 'ses-001', 'meg', ...
            'pmergedoe1000mspddfflo45hi45hfcstatic_001_array1.mat');
    end

    load(geomfile)
    D       = spm_eeg_load(datwithEMGmerged);
    grad_mm = D.sensors('MEG');
    ftdat   = spm2fieldtrip(D);

    badchans = D.chanlabels(D.badchannels);
    cfg = [];
    cfg.channel = setdiff(ftdat.label, badchans);
    ftdat = ft_selectdata(cfg, ftdat);

    if rectify
        cfg = [];
        cfg.rectify = 'yes';
        cfg.channel = 'EXG1';
        ftdatr = ft_preprocessing(cfg, ftdat);
        for k = 1:length(ftdat.trial)
            ftdat.trial{k}(end,:) = ftdatr.trial{k};
        end
    end

    %% load BEM v2 leadfield + position-based channel mapping
    load('C:\Users\mspedden\Documents\bem_spine_fields\bem_v2_leadfield_cervical_realistic_bem_.mat')

    nsourcepoints = size(leadfield_cord.pos, 1);
    spchanidx     = find(grad_mm.coilpos(:,2) < 200);
    spchanlabs    = grad_mm.label(spchanidx);

    % position-based mapping (v2 in metres, grad_mm in mm)
    grad_v2   = leadfield_cord.cfg.grad;
    spc_in_v2 = zeros(numel(spchanidx), 1);
    for i = 1:numel(spchanidx)
        p = grad_mm.coilpos(spchanidx(i), :);
        d = sqrt(sum((grad_v2.coilpos*1000 - p).^2, 2));
        [~, spc_in_v2(i)] = min(d);
    end

    Lf       = leadfield_cord;
    Lf.label = grad_mm.label(spchanidx);
    for i = 1:numel(leadfield_cord.leadfield)
        if ~isempty(leadfield_cord.leadfield{i})
            Lf.leadfield{i} = leadfield_cord.leadfield{i}(spc_in_v2, :);
        end
    end

    %% dummy head model
    cfg              = [];
    cfg.method       = 'infinite';
    cfg.siunits      = 1;
    cfg.grad         = grad_mm;
    cfg.conductivity = 1;
    dummyvol = ft_prepare_headmodel(cfg, mesh_torso);

    %% plot ROI on mesh (first subject only)
    if ss == 1
        figure('Color', 'w');
        ft_plot_mesh(mesh_wm, 'facecolor', [0.7 0.7 0.7], 'facealpha', 0.3, 'edgecolor', 'none');
        hold on
        highlightColor = [1 0.4 0.1];
        plot3(ROIpos(:,1), ROIpos(:,2), ROIpos(:,3), 'o', ...
            'MarkerSize', 10, ...
            'MarkerEdgeColor', [0.9 0.3 0], ...
            'MarkerFaceColor', highlightColor, ...
            'LineWidth', 2);
        view(90, 18);
        material dull;
    end

    %% LCMV beamformer
    cfg            = [];
    cfg.covariance = 'yes';
    tlock = ft_timelockanalysis(cfg, ftdat);

    cfg                       = [];
    cfg.method                = 'lcmv';
    cfg.headmodel             = dummyvol;
    cfg.sourcemodel.leadfield = Lf;
    cfg.grid                  = sources_cent;
    cfg.unit                  = Lf.unit;
    cfg.lcmv.keepfilter       = 'yes';
    source_idx = ft_sourceanalysis(cfg, tlock);

    %% virtual channel from ROI
    roi_center = mean(ROIpos, 1);
    d          = sqrt(sum((ROIpos - roi_center).^2, 2));
    R          = max(d);

    cfg              = [];
    cfg.pos          = roi_center;
    cfg.radius       = R;
    cfg.method       = 'svd';
    cfg.numcomponent = 1;
    VE = ft_virtualchannel(cfg, ftdat, source_idx);

    %% save without overwriting smoothed/unsmoothed versions
    savename = sprintf('VE_spine_sub%s_forspectra_bemv2%s', sub, ve_suffix);
    save(fullfile(save_dir, savename), 'VE')

end