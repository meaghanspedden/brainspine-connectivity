%Spinal cord virtual electrodes for FREQNESS


%% note will need to load both conditions and create filter based on both!

clear all
close all
clc

addpath('C:\Users\mspedden\Documents\brainspineconnectivity\source')
addpath('C:\Users\mspedden\Documents\spm')
spm('defaults','EEG')

addpath('C:\Users\mspedden\Documents\fieldtrip')
ft_defaults

save_dir='C:\Users\mspedden\Documents\brainspine_save';

%n=9
%subs = {'OP00212','OP00213',  'OP00215', 'OP00219', ...
   % 'OP00220', 'OP00221', 'OP00224', 'OP00225', 'OP00226'};

   sub='OP00212';

generic_dir = 'C:\Users\mspedden\Documents\new_leadfields_and_geom'; %where I have saved folder with brainspione leadfields and geoms (meshes)
geomfile = fullfile(generic_dir, 'geometries_cervical_realistic.mat');

LFop='spine'; %only want leadfields from spine here.


% subjResults=struct();
% 
% for ss=1:length(subs)
% 
%     sub=subs{ss};

    if strcmp(sub, 'OP00224') %saved under 002
%         dat = fullfile('C:\Users\mspedden\Documents', ...
%             ['sub-' sub], ...
%             'ses-001', ...
%             'meg', ...
%             'pmergedoe1000mspddfflo45hi45hfcstatic_002_array1.mat');

    else

        datcont = fullfile('C:\Users\mspedden\Documents', ...
            ['sub-' sub], ...
            'ses-001', ...
            'meg', ...
            'ddfflo45hi45hfcstatic_001_array1.mat');

        datrest= fullfile('C:\Users\mspedden\Documents', ...
            ['sub-' sub], ...
            'ses-001', ...
            'meg', ...
            'ddfflo45hi45hfcrest_002_array1.mat');
    end

    load(geomfile)
    Dcont=spm_eeg_load(datcont);
    Drest=spm_eeg_load(datrest);
    grad_mm=Dcont.sensors('MEG');
    ftdatcont = spm2fieldtrip(Dcont);
    ftdatrest = spm2fieldtrip(Drest);

    badchans=Dcont.chanlabels(Dcont.badchannels);

    %remove bad channels here.
    cfg=[];
    cfg.channel=setdiff(ftdatcont.label,badchans);
    ftdatcont=ft_selectdata(cfg,ftdatcont);
    ftdatrest=ft_selectdata(cfg,ftdatrest);


idx = startsWith(ftdatcont.label, {'X','Y','Z'});
OPMchans = ftdatcont.label(idx);

trllength=min([length(ftdatcont.trial{1}) length(ftdatrest.trial{1})]);

cfg=[];
cfg.begsample=1;
cfg.endsample=trllength;
ftdatcont=ft_redefinetrial(cfg,ftdatcont);
ftdatrest=ft_redefinetrial(cfg,ftdatrest);



ftdat=ft_appenddata([],ftdatcont,ftdatrest);

cfg                   = [];
cfg.covariance        = 'yes';
cfg.channel           = OPMchans;
tlock = ft_timelockanalysis(cfg, ftdat);

%% load and organise spinal cord leadfields
    [Gx, Gy, Gz] = build_leadfield_matrices(fullfile(generic_dir,'cervical_realistic_brain_spine'), LFop);

    nsourcepoints = size(Gx,1);
    nchannels     = size(Gx,2);
    spchanidx=find(grad_mm.coilpos(:,2) < 200); %indexed locally in grad struct (same as LF)
    spchanlabs=grad_mm.label(spchanidx);

    %% clip leadfields to spinal cord channels only
    Gx=Gx(:,spchanidx);
    Gy=Gy(:,spchanidx);
    Gz=Gz(:,spchanidx);
    
    %put leadfields into fieldtrip format
    Lf.pos    = sources_cent.pos;     % nsourcepoints x 3
    Lf.inside = sources_cent.inside;     % all points inside
    Lf.unit   = 'mm';
    Lf.label  = grad_mm.label(spchanidx);   % nchannels x 1 cell
    Lf.leadfielddimord = '{pos}_chan_ori';
    Lf.leadfield = cell(1,nsourcepoints);

    for k = 1:nsourcepoints
        % Combine X/Y/Z components like FT is used to
        Lf.leadfield{k} = [Gx(k,:)' Gy(k,:)' Gz(k,:)']; % nchannels x 3
    end

    % 2. dummy head model for input config only (not used)
    cfg                     = [];
    cfg.method              = 'infinite';
    cfg.siunits=1;
    cfg.grad=grad_mm;
    cfg.conductivity = 1;

    dummyvol = ft_prepare_headmodel(cfg,mesh_torso);

cfg                     = [];
cfg.method              = 'lcmv';
cfg.headmodel           = dummyvol;
cfg.sourcemodel.LF      = Lf.leadfield;
cfg.grid                = sources_cent;
cfg.unit                = Lf.unit;
cfg.lcmv.keepfilter     = 'yes';
cfg.channel             = OPMchans;
source_time = ft_sourceanalysis(cfg, tlock);

beamformer=source_time.avg.filter;

spine_source_data_cont = [];
spine_source_data_cont.label = arrayfun(@(x) sprintf('source%d', x), 1:nsourcepoints, 'UniformOutput', false);
spine_source_data_cont.time = ftdat.time;

spine_source_data_rest =[];
spine_source_data_rest.label =spine_source_data_cont.label;
spine_source_data_rest.time=ftdat.time;


for sp=1:nsourcepoints
    for i=1:length(ftdatcont.trial) %gives a time series for each ori
        spine_source_data_cont.trial{i} = beamformer{sp} * ftdatcont.trial{i}(idx,:);
        spine_source_data_rest.trial{i} = beamformer{sp} * ftdatrest.trial{i}(idx,:);
    end
end


timeseries_c = cat(2, spine_source_data_cont.trial{:});
[u1, s1, v1] = svd(timeseries_c, 'econ');

timeseries_r = cat(2,spine_source_data_cont.trial{:});
[u2, s2, v2] = svd(timeseries_r, 'econ');


virtualchanneldata_c = [];
virtualchanneldata_c.label = spine_source_data_cont.label;
virtualchanneldata_c.time = ftdatcont.time;

virtualchanneldata_r = [];
virtualchanneldata_r.label = spine_source_data_cont.label;
virtualchanneldata_r.time = ftdatcont.time;


    for k = 1:length(ftdatcont.trial)
        for sp=1:nsourcepoints
        virtualchanneldata_c.trial{k}(sp,:) = u1(:,1)' * beamformer{sp} * ftdatcont.trial{k}(idx,:);
        virtualchanneldata_r.trial{k}(sp,:) = u2(:,1)' * beamformer{sp} * ftdatrest.trial{k}(idx,:);

        end
    end


savename=sprintf('datVC_%s',sub);
datmat1=virtualchanneldata_c.trial{1};
datmat2=virtualchanneldata_r.trial{1};
save(fullfile('C:\Users\mspedden\Documents\Frequency-resolved_brain_network_estimation_via_source_separation_FREQ-NESS\FREQNESS_Toolbox\FREQNESS_Data\Dataset_1',savename),'datmat1')
save(fullfile('C:\Users\mspedden\Documents\Frequency-resolved_brain_network_estimation_via_source_separation_FREQ-NESS\FREQNESS_Toolbox\FREQNESS_Data\Dataset_2',savename),'datmat2')
