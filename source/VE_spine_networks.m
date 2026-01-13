% virtual electrode spinal cord sources (all) 
clear all
close all
clc

addpath('C:\Users\mspaedden\Documents\brainspineconnectivity\source')
addpath('C:\Users\mspedden\Documents\spm')
spm('defaults','EEG')

addpath('C:\Users\mspedden\Documents\fieldtrip')
ft_defaults

save_dir='C:\Users\mspedden\Documents\brainspine_save';
if ~exist(save_dir,'dir')
    mkdir(save_dir)
end


subs = {'OP00212'};%,'OP00213','OP00215', 'OP00219', ...
% 'OP00225', 'OP00221', 'OP00224'}; %have brain 

generic_dir = 'C:\Users\mspedden\Documents\new_leadfields_and_geom'; %where I have saved folder with brainspione leadfields and geoms (meshes)
geomfile = fullfile(generic_dir, 'geometries_cervical_realistic.mat');

LFop='spine';
HFC=1;
fband=[5 45];



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

    [Gx, Gy, Gz] = build_leadfield_matrices(fullfile(generic_dir,'cervical_realistic_brain_spine'), LFop);

    nsourcepoints = size(Gx,1);
    nchannels     = size(Gx,2);
    spchanidx=find(grad_mm.coilpos(:,2) < 200); %indexed in grad struct
    spchanlabs=grad_mm.label(spchanidx);


 %% clip leadfields to spinal cord channels only
    Gx=Gx(:,spchanidx);
    Gy=Gy(:,spchanidx);
    Gz=Gz(:,spchanidx);

    Lf.pos    = sources_cent.pos;     % nsourcepoints x 3
    Lf.inside = sources_cent.inside;     % all points inside
    Lf.unit   = 'mm';
    Lf.label  = grad_mm.label(spchanidx);   % nchannels x 1 cell
    Lf.leadfielddimord = '{pos}_chan_ori';

    Lf.leadfield = cell(1,nsourcepoints);

    for k = 1:nsourcepoints
        % Combine X/Y/Z components
        Lf.leadfield{k} = [Gx(k,:)' Gy(k,:)' Gz(k,:)']; % nchannels x 3
    end

    % 2. dummy head model for input config only
    cfg                     = [];
    cfg.method              = 'infinite';
    cfg.siunits=1;
    cfg.grad=grad_mm;
    cfg.conductivity = 1;

    dummyvol = ft_prepare_headmodel(cfg,mesh_torso);
   

%% lcmv

cfg                   = [];
cfg.covariance        = 'yes';
tlock = ft_timelockanalysis(cfg, ftdat);

cfg                     = [];
cfg.method              = 'lcmv';
cfg.headmodel           = dummyvol;
cfg.sourcemodel.leadfield     = Lf;
cfg.grid     = sources_cent;
cfg.unit                = Lf.unit;
cfg.lcmv.keepfilter     = 'yes';
source_idx = ft_sourceanalysis(cfg, tlock);

cfg=[];
cfg.pos=sources_cent.pos;
cfg.method='pca';
VE = ft_virtualchannel(cfg, ftdat, source_idx);

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

savename=sprintf('VE_spine_sub%s_allsources', sub);
save(fullfile(save_dir,savename), 'VE_s', 'VE_r')


end