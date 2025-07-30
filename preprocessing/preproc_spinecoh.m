%% process data from hand contraction. ~ July 2025


clear all;
close all;

sub='OP00212';

rmchans=0; %remove bad channels from psd...for first round of analysis

analysis={'static','rest'};
save_dir = fullfile('D:\MSST001', [sub '_merged']);

if ~exist(save_dir, 'dir')
    mkdir(save_dir);
end

%% heartbeat estimate from 'best segment'
if strcmp(sub, 'OP00220')

    flag=1;
else
    flag=0;
end

%% paths (need to tidy this)
addpath('D:\spm')
spm('defaults','EEG')
addpath('D:\brainspineconnectivity\plotting')
addpath('D:\brainspineconnectivity\preprocessing\version_2')
addpath('D:\msg_source_recon')
addpath('D:\brainspineconnectivity')

%% labels
brainchan_labels={'H1', 'H6', 'B3', 'H8', 'G3', 'B6', 'H7', 'H4', 'G4', 'B4', 'B2', 'H2', 'B7', 'B8', 'B5', 'B1'};
heart_labels={'C1', 'D3', 'D2', 'D4'};

spine_labels={'C2', 'C3', 'C4', 'C7', 'D8', 'D5', 'D6', 'C5', 'C8', 'C6', 'D7', ...
    'F2', 'F7', 'F5', 'F8', 'F3', 'F1', 'F4', 'F6', 'E2', 'E4', 'E8', 'E5', ...
    'E1', 'E7', 'E6', 'E3', 'H5', 'G5', 'G7', 'H3', 'G1', 'G2', 'A7', 'G6', ...
    'G8', 'A6', 'A3', 'A2', 'A8', 'A5', 'A1', 'A4'};
trigger_label='T5';

filetemplate = fullfile('D:\MSST001', ['sub-' sub], 'ses-001', 'meg');

posfile=fullfile(filetemplate, 'static_001_ar_positions.tsv');

EMGpath = fullfile('D:\', [sub '_experiment'], 'EMG');

EMGfiletemplate = [sub '_00'];


%% Processing Options
fband=[5 45]; %% filter band
stpband1=[48 52];
stpband2=[98 102];
EpochTime=1; %  %% length of epoch in seconds
DataTime=115; %length to crop opm and emg to, s

% flags to turn these filters on/off
HFCflag=0; %% homogenous field correction (0 or 1)
HB=1; %% heartbeat removal


%% preproc loop

allfilenames={};

for cond=1:length(analysis)

    [runs, ref_labels, heartidx, badchanlabels] = get_metadata(sub,analysis{cond}); %only runs differ here woul be better to output across condition


    for k=1:numel(runs.opm)

        S=[];
        S.data = fullfile(filetemplate, sprintf('%s_%s_array1.lvm', analysis{cond},runs.opm{k}));
        S.positions = posfile;
        D = spm_opm_create(S);


        %% psd
        opms=D.indchannel(D.sensors('MEG').label);

        S = [];
        S.D = D;
        S.plot = 1;
        S.channels = D.chanlabels(opms);
        S.triallength = 2000;
        S.wind = @hanning;
        
        if rmchans && k==1 && cond==1
            S.selectbad=1;
            [~,~,badidx] = spm_opm_psd(S);

            if ~isempty(ref_labels)
                refidx=find(contains(D.chanlabels,ref_labels));
                badidx=[badidx refidx]; %for now set ref as bad
            end

            if~isempty(badchanlabels)
                badsens=find(contains(D.chanlabels,badchanlabels));
                badidx=[badidx badsens];
            end

            save(fullfile(save_dir,sprintf('%s_badchans',sub)), 'badidx')

        else 
            S.selectbad=0; 
            spm_opm_psd(S);
            load(fullfile(save_dir,sprintf('%s_badchans',sub)), 'badidx') %same bad channels across runs
       end

        D=badchannels(D, badidx,1);


        %% read in EMG data

        filename=sprintf('%s_00%s.vhdr', sub, runs.emg{k});

        cfg = [];
        cfg.dataset = fullfile(EMGpath,filename);

        EMGdata = ft_preprocessing(cfg);

        EMG_trigs=ft_read_event(fullfile(EMGpath,sprintf('%s_00%s.eeg', sub, runs.emg{k})));

        D_EMG = spm_eeg_ft2spm(EMGdata, sprintf('EMG_%s', sub));


        if HFCflag
            S=[];
            S.prefix='hfc';
            S.D=D;
            fprintf('\n Running HFC');
            [D_hfc, ~]=spm_opm_hfc(S);

            S = [];
            S.D1 =D_hfc;
            S.D2=D;
            S.plot = 1;
            %S.channels = D.chanlabels(setdiff(opms,badidx));
            S.triallength = 2000;
            S.wind = @hanning;
            [shield,~] = spm_opm_rpsd(S);

            D=D_hfc;
        end
        %%  Now filter

        S=[];
        S.D=D;
        S.freq=min(fband);
        S.order=5;
        S.band='high';
        S.prefix=sprintf('hi%d',fband(2));
        D=spm_eeg_filter(S);

        S.D=D_EMG;
        D_EMG=spm_eeg_filter(S);

        S=[];
        S.D=D;
        S.freq=max(fband);
        S.order=5;
        S.band='low';
        S.prefix=sprintf('lo%d',fband(2));
        D=spm_eeg_filter(S);

        S.D=D_EMG;
        D_EMG=spm_eeg_filter(S);


        S=[];
        S.D=D;
        S.band='stop';
        S.order=3;
        S.freq=stpband1;
        S.D=D;
        D=spm_eeg_filter(S);

        S.D=D_EMG;
        D_EMG=spm_eeg_filter(S);

        S=[];
        S.D=D;
        S.band='stop';
        S.order=3;
        S.freq=stpband2;
        D=spm_eeg_filter(S);

        S.D=D_EMG;
        D_EMG=spm_eeg_filter(S);

        heart_labels_1 = labelToIndex(heart_labels);
        heartind=find(contains(D.chanlabels, heart_labels_1));

        megind=setdiff(D.indchantype('MEG'), D.badchannels);


        if HB % estimate heartbeat over all channels

            if k==1 && cond==1
                [heartest,beatlen]=est_heartbeat(D,heartind,megind,heartidx,flag);
                allheart=zeros(length(cond), numel(runs.opm),length(megind),beatlen);
            else
                [heartest]=est_heartbeat(D,heartind,megind,heartidx,flag,beatlen);
            end
            allheart(cond, k, :, :) = heartest;
        end % if HB


        %% downsample OPM data to match EMG

        S=[];
        S.fsample_new=fsample(D_EMG);
        S.D=D;
        opm_crop=spm_eeg_downsample(S);

if ~isempty(ref_labels)
        refidx=find(contains(opm_crop.chanlabels,ref_labels));
        opm_crop=chantype(opm_crop,refidx,'REF');

        S=[];
        S.D=opm_crop;
        S.confounds={'REF'};

        opm_crop = spm_opm_synth_gradiometer(S);
end

        % find peak for trigger opm

        trigger_idx=find(contains(D.chanlabels, trigger_label));
        %find when trig comes
        eventCont           = D(trigger_idx,:);
        time                = D.time;
        sample              = 1:length(time);
        % convert to binary
        eventDisc           = (eventCont > 0.9);

        % find first sample. first find differences
        tmp                 = eventDisc - [0,eventDisc(1:end-1)];
        tmp2                = (tmp == 1);
        events              = sample(tmp2)';
        events              = round(events);

        %because it is downsampled
        event = round((events/(D.fsample/opm_crop.fsample)));


        %crop opms
        S = [];
        S.D = opm_crop;
        S.timewin=[opm_crop.time(event(1))*1000 (opm_crop.time(event(1))+DataTime)*1000]; %ms!!!
        opm_crop = spm_eeg_crop(S);

        %crop emg
        EMG_start=D_EMG.time(EMG_trigs(2).sample); %first trigger in this system is a start trigger

        S=[];
        S.D=D_EMG;
        S.timewin=[EMG_start*1000 (EMG_start+DataTime)*1000];
        emg_crop=spm_eeg_crop(S);

        EMGchanname='EXG1';
        chan2replace='Data_Valid2';

        %add emg to opm object
        replaceidx=opm_crop.indchannel(chan2replace);

        opm_crop(replaceidx,:)=opm_crop(replaceidx,:).*0;
        opm_crop(replaceidx,:)=emg_crop(1,:);
        opm_crop = chanlabels(opm_crop,replaceidx,EMGchanname);
        opm_crop = units(opm_crop,replaceidx,'uV');
        opm_crop = chantype(opm_crop,replaceidx,'EMG');

        opm_crop.save;


        %epoch into 1-s trials

        epochsamples=EpochTime*opm_crop.fsample;
        nsamples=size(opm_crop(:,:),2);

        Nepochs=floor(nsamples/epochsamples);

        %make trl matrix for epoching
        trl=zeros(Nepochs,3);
        for f=1:Nepochs
            trl(f,:)=[(f-1)*epochsamples+1 f*epochsamples 0];
        end


        condition_labels = cell(Nepochs,1);
        condition_labels(:) = {analysis{cond}};

        S=[];
        S.D=opm_crop;
        S.trl=trl;
        S.bc=0;
        S.conditionlabels=condition_labels;
        S.prefix=sprintf('e%dms',round(EpochTime*1000));
        Dep = spm_eeg_epochs(S);


        %% remove outliers
        S=[];
        S.D=Dep;
        [Dep, retain] = spm_opm_removeOutlierTrials(S);

        savename=sprintf('retainedtrials%srun%s_%s',sub,runs.opm{k},analysis{cond});
        save(fullfile(save_dir,savename),'retain')

        allfilenames=strvcat(allfilenames,Dep.fname);

        ftdat=spm2fieldtrip(Dep);
        ftdat=rmfield(ftdat, 'hdr');
        badidx=D.badchannels;
        badlabs=ftdat.label(badidx);
        labels = ftdat.label;
        isGood = ~contains(labels, badlabs);
        hasXYZ = contains(labels, 'X') | contains(labels, 'Y') | contains(labels, 'Z');
        goodlabs = labels(isGood & hasXYZ);
        cfg=[];
        cfg.channel=goodlabs;
        cfg.allowoverlap='yes';
        ft_databrowser(cfg,ftdat)


    end % for runs
end %loop through conditions
% now merge (before heartbeat removal)


cd(opm_crop.path);
S=[];
S.D=allfilenames;
S.prefix='merged';
DAll=spm_eeg_merge(S);

if HB==1
    avheart = squeeze(mean(allheart, [1 2]));  % mean over cond and k
    figure
    plot(avheart')

    BALANCE=1;
    Dbeforeproject=DAll;
    DAll=remove_heartbeat(DAll,avheart,megind,BALANCE);
end



% [coh] = getCoh(DAll,brainchan_labels);


%final merged file

