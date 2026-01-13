% virtual electrode brain - network stuff
clear all
close all
clc

addpath('C:\Users\mspedden\Documents\brainspineconnectivity\source')
addpath('C:\Users\mspedden\Documents\spm')
spm('defaults','EEG')

addpath('C:\Users\mspedden\Documents\fieldtrip')
ft_defaults

save_dir='C:\Users\mspedden\Documents\brainspine_save';
if ~exist(save_dir,'dir')
    mkdir(save_dir)
end

subs={'OP00212'};

% subs = {'OP00212','OP00213','OP00215', 'OP00219', ...
%  'OP00225', 'OP00221', 'OP00224'}; %have brain 

generic_dir = 'C:\Users\mspedden\Documents\new_leadfields_and_geom'; %where I have saved folder with brainspione leadfields and geoms (meshes)
geomfile = fullfile(generic_dir, 'geometries_cervical_realistic.mat');

HFC=1;

for ss=1:length(subs)

    sub=subs{ss};

%     if strcmp(sub, 'OP00224')
%         datwithEMGmerged = fullfile('C:\Users\mspedden\Documents', ...
%             ['sub-' sub], ...
%             'ses-001', ...
%             'meg', ...
%             'pmergedoe1000mspddfflo45hi45hfcstatic_002_array1.mat');
%     else
        datrest = fullfile('C:\Users\mspedden\Documents', ...
            ['sub-' sub], ...
            'ses-001', ...
            'meg', ...
            'ddfflo45hi45hfcrest_002_array1.mat');

        datstat = fullfile('C:\Users\mspedden\Documents', ...
            ['sub-' sub], ...
            'ses-001', ...
            'meg', ...
            'ddfflo45hi45hfcstatic_001_array1.mat');

    %end


    load(geomfile)

    Dr=spm_eeg_load(datrest);
    Ds=spm_eeg_load(datstat);
    grad_mm=Dr.sensors('MEG');
    ftdat_r = spm2fieldtrip(Dr);
    ftdat_s=spm2fieldtrip(Ds);

    len_r=length(ftdat_r.trial{1}); %to separate later
    len_s=length(ftdat_s.trial{1});
    
    %concatenate conditions for constructing bf filter
    ftdat = ftdat_r;

    fs = ftdat.hdr.Fs;

    ftdat.trial{1} = [ftdat_r.trial{1}  ftdat_s.trial{1}];
    ftdat.time{1}  = [ftdat_r.time{1} ...
        ftdat_s.time{1} + ftdat_r.time{1}(end) + 1/fs];

    ftdat.trialinfo = [];
   
    badchans=Ds.chanlabels(Ds.badchannels);

    %remove bad channels
    cfg=[];
    cfg.channel=setdiff(ftdat.label,badchans);
    ftdat=ft_selectdata(cfg,ftdat);

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

%figure; ft_plot_sens(braingrad)

braindat=ftdat;

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
LF = ft_prepare_leadfield(cfg);


%% lcmv

cfg                   = [];
cfg.covariance        = 'yes';
tlock = ft_timelockanalysis(cfg, braindat);

cfg                     = [];
cfg.method              = 'lcmv';
cfg.headmodel           = vol;
cfg.sourcemodel.LF     = LF;
cfg.grid     = sources_brain;
cfg.unit                = LF.unit;
cfg.lcmv.keepfilter     = 'yes';
source_idx = ft_sourceanalysis(cfg, tlock);

cfg=[];
cfg.pos=source_idx.pos;
cfg.method='pca'; %dont think I need this
VE = ft_virtualchannel(cfg, braindat, source_idx);

%split back into two conditions
VE_r = VE;
VE_s = VE;

idx_r = 1:len_r;
idx_s = len_r+1 : size(VE.trial{1},2);

% split data
VE_r.trial{1} = VE.trial{1}(:, idx_r);
VE_s.trial{1} = VE.trial{1}(:, idx_s);

VE_r.time{1} = VE.time{1}(idx_r);
VE_s.time{1} = VE.time{1}(idx_s);

VE_r.trialinfo = [];
VE_s.trialinfo = [];

savename=sprintf('VE_brain_sub%s_allsources', sub);
save(fullfile(save_dir,savename), 'VE_s', 'VE_r')


end