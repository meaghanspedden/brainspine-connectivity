
%% source recon fieldtrip brain

clear all
close all
clc

addpath('C:\Users\mspedden\Documents\brainspineconnectivity\source')
addpath('C:\Users\mspedden\Documents\spm')
spm('defaults','EEG')

addpath('C:\Users\mspedden\Documents\fieldtrip')
ft_defaults

subs={'OP00212'};

%subs = {'OP00212','OP00213',  'OP00215', 'OP00219', ...
% 'OP00220', 'OP00221', 'OP00224', 'OP00225', 'OP00226'};

generic_dir = 'C:\Users\mspedden\Documents\new_leadfields_and_geom'; %where I have saved folder with brainspione leadfields and geoms (meshes)
geomfile = fullfile(generic_dir, 'geometries_cervical_realistic.mat');

HFC=1;
%LFop='spine';
correct_multcomp=1;

subjResults=struct();

for ss=1:length(subs)

    sub=subs{ss};

    if strcmp(sub, 'OP00224')
        datwithEMGmerged = fullfile('C:\Users\mspedden\Documents', ...
            ['sub-' sub], ...
            'ses-001', ...
            'meg', ...
            'pmergedoe1000mspddfflo45hi45hfcstatic_002_array1.mat');
    else
        datwithEMGmerged = fullfile('C:\Users\mspedden\Documents', ...
            ['sub-' sub], ...
            'ses-001', ...
            'meg', ...
            'pmergedoe1000mspddfflo45hi45hfcstatic_001_array1.mat');
    end

    load(geomfile)
    D=spm_eeg_load(datwithEMGmerged);
    grad_mm=D.sensors('MEG');
    ftdat = spm2fieldtrip(D);

    badchans=D.chanlabels(D.badchannels);

    %remove bad channels here.
    cfg=[];
    cfg.channel=setdiff(ftdat.label,badchans);
    ftdat=ft_selectdata(cfg,ftdat);

    %% try simple volume conductor first


    %     [Gx, Gy, Gz] = build_leadfield_matrices(fullfile(generic_dir,'cervical_realistic_brain_spine'), LFop);
    %
    %     nsourcepoints = size(Gx,1);
    %     nchannels     = size(Gx,2);

    %brainchanidx=find(grad_mm.coilpos(:,2) > 200);
    
    brainlabs=[];
    for f = 1:length(grad_mm.label)
        if grad_mm.chanpos(f,2) > 200
            brainlabs = [brainlabs grad_mm.label(f)];
        end
    end

    %% new brain only dataset and grad

brainidx = find(grad_mm.chanpos(:,2) > 200);
braingrad = grad_mm;
braingrad.chanori = grad_mm.chanori(brainidx, :);
braingrad.chanpos = grad_mm.chanpos(brainidx, :);
braingrad.chantype = grad_mm.chantype(1:length(brainidx));
braingrad.chanunit = grad_mm.chanunit(1:length(brainidx));
braingrad.coilori = grad_mm.coilori(brainidx, :);
braingrad.coilpos = grad_mm.coilpos(brainidx, :);
braingrad.label = grad_mm.label(brainidx);
braingrad.tra = grad_mm.tra(brainidx, brainidx);

figure; ft_plot_sens(braingrad)


cfg=[];
cfg.channel=brainlabs; %only use brain channels
braindat=ft_selectdata(cfg,ftdat);

mesh_brain.unit='mm';
%single shell volume conductor
cfg = [];
cfg.method = 'singleshell';
vol = ft_prepare_headmodel(cfg, mesh_brain);

% calculate forward
cfg                     = [];
cfg.sourcemodel         = sources_brain;
cfg.headmodel           = vol;
cfg.grad                = braingrad;
cfg.reducerank          = 'no';
%cfg.normalize           = 'yes';
%cfg.normalizeparam     = 0.5;

LF = ft_prepare_leadfield(cfg);


    %% beamforming----------------------------
    %1. trial wise freq dat
    cfg=[];
    cfg.output     = 'powandcsd';
    cfg.method     = 'mtmfft';
    cfg.foilim     = [10 35];
    cfg.tapsmofrq  = 1;
    cfg.keeptrials='yes';
    freqdat_tr=ft_freqanalysis(cfg,braindat); %trial wise freq

    %% separate conditions

    statidx=find(ftdat.trialinfo==1);
    restidx=find(ftdat.trialinfo==2);

    cfg=[];
    cfg.trials=statidx;
    statdat=ft_selectdata(cfg,freqdat_tr); %trial wise data pr condition

    cfg.trials=restidx;
    restdat=ft_selectdata(cfg,freqdat_tr);

    %% get struct with mean freq data for constructing bf filter
    cfg=[];
    cfg.output     = 'powandcsd';
    cfg.method     = 'mtmfft';
    cfg.foilim     = [10 35];
    cfg.tapsmofrq  = 1;
    cfg.keeptrials = 'no';
    freqdat=ft_freqanalysis(cfg,braindat);

    %% construct filter using both conditions, use avg data
    cfg=[];
    cfg.grid = sources_brain;
    cfg.headmodel=vol;
    cfg.sourcemodel.leadfield=LF;
    cfg.dics.keepfilter='yes';
    cfg.dics.lambda=10;
    cfg.method = 'dics';
    %cfg.fixedori='yes';
    source_all = ft_sourceanalysis(cfg,freqdat);

    %% get trial-wise estimates with filter from above per condition
    cfg=[];
    cfg.grid = sources_brain;
    cfg.headmodel=vol;
    cfg.sourcemodel.leadfield=LF;
    cfg.dics.filter=source_all.avg.filter;
    cfg.dics.lambda=10;
    cfg.method = 'dics';
    cfg.rawtrial='yes';
    cfg.keeptrials='yes';

    source_stat = ft_sourceanalysis(cfg,statdat); %struct per condition with single trials
    source_rest = ft_sourceanalysis(cfg,restdat);

%% --- Convert trial structs to matrices ---
nSources = numel(source_stat.inside);

nTrials_task = numel(source_stat.trial);
pow_task = nan(nSources, nTrials_task);
for t = 1:nTrials_task
    pow_task(:,t) = source_stat.trial(t).pow;
end

nTrials_rest = numel(source_rest.trial);
pow_rest = nan(nSources, nTrials_rest);
for t = 1:nTrials_rest
    pow_rest(:,t) = source_rest.trial(t).pow;
end

nTrials_min = min(nTrials_task, nTrials_rest);
pow_task = pow_task(:,1:nTrials_min);
pow_rest = pow_rest(:,1:nTrials_min);

source_task_clean = [];
source_task_clean.pow    = pow_task;
source_task_clean.inside = source_stat.inside;
source_task_clean.pos    = source_stat.pos;
source_task_clean.unit   = source_stat.unit;

source_rest_clean = [];
source_rest_clean.pow    = pow_rest;
source_rest_clean.inside = source_rest.inside;
source_rest_clean.pos    = source_rest.pos;
source_rest_clean.unit   = source_rest.unit;

source_task_clean.pow = log10(source_task_clean.pow);
source_rest_clean.pow = log10(source_rest_clean.pow);
%% do stats
cfg = [];
%cfg.method    = 'analytic';       % simple t-test (or 'montecarlo' for cluster-based)
cfg.statistic = 'depsamplesT';    % paired t-test
cfg.parameter = 'pow';
%cfg.correctm  = 'no';             % no multiple comparison correction

cfg.method            = 'montecarlo';           % permutation-based test
cfg.statistic         = 'depsamplesT';          % still paired
cfg.correctm          = 'cluster';             % or 'no' for uncorrected
cfg.clusteralpha      = 0.05;                   % threshold for forming clusters
cfg.clusterstatistic  = 'maxsum';               % cluster statistic
cfg.tail              = 0;                      % two-tailed
cfg.clustertail       = 0;                      % two-tailed clusters
cfg.alpha             = 0.05;                   % significance level
cfg.numrandomization  = 1000;                   % number of permutations



% Design: row1 = subject index (same for all trials, single subject)
%         row2 = condition label (1=task, 2=rest)
design = [ones(1,nTrials_min*2); [ones(1,nTrials_min), 2*ones(1,nTrials_min)]];
cfg.design = design;
cfg.uvar   = 1; 
cfg.ivar   = 2; 

stat = ft_sourcestatistics(cfg, source_task_clean, source_rest_clean);

% %% interpolate on MRI

% mri = ft_read_mri('C:\Users\mspedden\Documents\new_leadfields_and_geom\padded_combined_image_1point5iso.nii'); 
% 
% cfg = [];
% cfg.parameter    = 'stat';
% cfg.maskparameter = 'mask';
% cfg.interpmethod = 'nearest';
% source_int = ft_sourceinterpolate(cfg, stat_clean, mri);
% % 
% cfg_plot = [];
% cfg_plot.method       = 'ortho';     % orthogonal slices
% cfg_plot.funparameter = 'stat';
% cfg_plot.maskparameter= 'stat';       % uses NaNs from interpolation
% cfg_plot.funcolormap  = 'jet';
% %cfg_plot.funcolorlim  = [-5 5];      % adjust t-value range
% ft_sourceplot(cfg_plot, source_int);



%% interpolate on brain mesh
brainmesh=struct();
brainmesh.pos=mesh_brain.vertices;
brainmesh.tri=mesh_brain.faces;

stat_clean = [];
stat_clean.stat   = stat.stat;
stat_clean.mask   = stat.mask;
stat_clean.pos    = stat.pos;
stat_clean.inside = stat.inside;
stat_clean.unit   = stat.unit;  % optional

cfg = [];
cfg.parameter    = 'stat';        % t-values
cfg.maskparameter = 'mask';        % only interpolate significant points
cfg.interpmethod = 'nearest';
source_int = ft_sourceinterpolate(cfg, stat_clean, brainmesh);

%%
cfg_plot = [];
cfg_plot.method        = 'surface';
cfg_plot.funparameter  = 'stat';        % masked t-values
cfg_plot.maskparameter = 'mask';        % use the NaNs from interpolation as a mask
cfg_plot.funcolormap   = 'jet';
%cfg_plot.funcolorlim   = [-2 2];       % adjust to your t-value range
%cfg_plot.opacitymap    = 'rampup';     % transparent non-significant areas
cfg_plot.projmethod    = 'nearest';    % projection onto surface
ft_sourceplot(cfg_plot, source_int);
view(179,-43)
camlight;
lighting gouraud;


end




