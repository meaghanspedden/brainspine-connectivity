%% Pre process data from hand contraction. ~ July 2023
% concurrent brain spinal cord OPM recordings

%subject 136

clear all;
close all;

sub='OP00136';

bids_session='001';
dsprefix=1;

brainchan_labels={'MN','OI', 'MR', 'A4','K5', 'OG','DU','DH','DL','DO','35','MF','N2','MZ','MU','MK'};
badchans={'KR-X', 'KR-Y', 'KR-Z', 'MJ-X', 'MJ-Y', 'MJ-Z'}; %%
posfileneck='D:\MSST001\sub-OP00136\ses-001\meg\ds_sub-OP00136_ses-001_positions.tsv';

EMGpath='D:\OP00136_experiment\EMGfiles';
EMGfiletemplate='000136_static_';
savedir='D:\MSST001\Coh_results00136';



%% Processing Options
    fband=[5 45]; %% OPM filter band
    stpband1=[48 52];
    stpband2=[98 102];
    emghpf=10; %% hp for emg
    EpochTime=1; %  %% length of epoch in seconds

    % flags to turn these filters on/off
    HFCflag=0; %% homogenous field correction (0 or 1)
    MEGHP=1;
    MEGLP=1;
    MEGBS=1; % stop band
    HB=1; %% heartbeat removal
    EMGHP=1; %% HP filter EMG
    ONEEPOCH=0; %% put this on if you just want single long epoch


%% which opm files correspond to which emg files
Bothexptorder={...
    '002'  'R'  '08'
    '003'  'R'  '10';
    '005'  'R'  '12';
    '006'  'R'  '14';

    '002'  'L'  '16';
    '003'  'L'  '18';
    '004'  'L'  '20';
    '005'  'L'  '22'};


%% trig info

OPMEMGsync='NI-TRIG-1'; %% EMG trigger syncing to OPM
OPMchansreplace=strvcat('NI-TRIG-2','NI-TRIG-3'); %% replace this/these OPM channels with EMG data
EMGtriglatency=0.05;

%% paths

datadir= 'D:\MSST001';
addpath('D:\spm')
spm('defaults','EEG')
addpath(genpath('D:\brainspineconnectivity'))
warning('MAKE SURE TO CHANGE spm_eeg_filter to GRB version')

%% right and left hands

for RIGHT=1

    if RIGHT
        exptorder=Bothexptorder(1:4,:);
        pstr='rhand';
        hbcomps=[2 3 3 3];
    else
        exptorder=Bothexptorder(5:8,:);
        pstr='lhand';
        hbcomps=[2 2 2 1]; %components containing heartbeat activity in each run. hard coded.
    end

    %%

    allfilenames=[];
    allxfilenames=[];

    for exptind=1:size(exptorder,1) %loop through runs

        if strcmp(exptorder(exptind,2),'R')
            EMGchanname='EMG 1';
            otherEMGchan='EMG 2';
            EMGfilename=[EMGfiletemplate,cell2mat(exptorder(exptind,3)),'.smrx'];
            task='staticright';
        end

        if strcmp(exptorder(exptind,2),'L')
            EMGchanname='EMG 2';
            otherEMGchan='EMG 1';
            EMGfilename=[EMGfiletemplate,cell2mat(exptorder(exptind,3)),'.smrx'];
            task='staticleft';
        end

        run=cell2mat(exptorder(exptind,1));

        %% make SPM object
         [D,fullfilename]=convert_opm2spm(datadir,posfileneck,sub,bids_session,task,run,'single',dsprefix);
        
         %for 136 1b and ds come in and out
        if strcmp(sub,'OP00136') 
            chanlabs=D.chanlabels;
            S=[];
            S.channels=chanlabs(~contains(chanlabs,{'G2-1B-Y','G2-1B-Z','G2-DS-Y','G2-DS-Z'}));
            S.D=D;
            D=spm_eeg_crop(S);
        end

        %% psd
        opms=D.indchannel(D.sensors('MEG').label);

        S = [];
        S.D = D;
        S.plot = 1;
        S.channels =D.chanlabels(opms);
        S.triallength = 2000;
        S.wind = @hanning;
        spm_opm_psd(S);
        xlim([1,100])
        title(sprintf('raw psd run %g',exptind))

        %% visual inspection time series for bad channels

        labs=D.chanlabels;
        datsamp=D(opms,1:4000);
        figure; plot(datsamp')
        legend(labs(opms))

        %% read in EMG data
        [~, EMGstruct]= readSpikeData(fullfile(EMGpath, EMGfilename),EMGchanname);

        cfg = [];
        cfg.output    = 'pow';
        cfg.channel   = 'all';
        cfg.method    = 'mtmfft';
        cfg.taper     = 'dpss';
        cfg.pad       ='maxperlen';
        cfg.tapsmofrq = 3;
        cfg.foi       = 0:1:100;
        emg_freq    = ft_freqanalysis(cfg, EMGstruct);

        figure; semilogy(emg_freq.freq, emg_freq.powspctrm)
        xlabel('Frequency'); ylabel('Log power')
        title('EMG power spectrum')


        if MEGHP
            S=[];
            S.D=D;
            S.freq=min(fband);
            S.order=5;
            S.band='high';
            S.chans2filter=D.indchantype('MEG');
            S.prefix=sprintf('hi%d',fband(1));
            D=spm_eeg_filter(S);


        datft=spm2fieldtrip(D);
        datft=rmfield(datft,'hdr');
        cfg=[];
        cfg.channel=datft.grad.label;
        %dbfig=ft_databrowser(cfg,datft)


%         S=[];
%         S.timewin=[0 135*1000]; %in ms %need to check this for all trials
%         S.D=D;
%         D=spm_eeg_crop(S);

        end
        %%  NOW LP FILTER JUST MEG
        if MEGLP
            S=[];
            S.D=D;
            S.freq=max(fband);
            S.order=5;
            S.band='low';
            S.prefix=sprintf('lo%d',fband(2));
            S.chans2filter=D.indchantype('MEG'); %% JUSTMEG NOT EMG
            D=spm_eeg_filter(S);
        end


        %% BAND STOP MEG CHANS
        if MEGBS
            S=[];
            S.D=D;
            S.band='stop';
            S.order=3;
            S.freq=stpband1;
            S.D=D;
            S.chans2filter=D.indchantype('MEG');
            D=spm_eeg_filter(S);

            S=[];
            S.D=D;
            S.band='stop';
            S.order=3;
            S.freq=stpband2;
            S.chans2filter=D.indchantype('MEG'); % don't apply 2nd bandstop to emg
            D=spm_eeg_filter(S);

        end

%         datft=spm2fieldtrip(D);
%         datft=rmfield(datft,'hdr');
%         cfg=[];
%         cfg.channel=datft.grad.label;
%         dbfig=ft_databrowser(cfg,datft)
%% get spine channel indices and channel labels

        spineind=D.indchannel(D.sensors('MEG').label);
        allspineind(exptind,:)=spineind;
        allchannames(exptind,:)=D.sensors('MEG').label;


        %% mark bad channels
        badind=[];
        for f=1:length(badchans)
           badidx= find(contains(D.chanlabels,badchans{f}));
           if ~isempty(badidx)
            badind(f)=badidx;
           end
        end
        D = badchannels(D, badind, 1); %% set channels to bad

    

        %% heartbeat estimation
        if HB   %% estimate heartbeat over all channels (will remove after merging files)
             
             if exptind==1
                [heartest,beatlen,megind]=grb_est_heartbeat(D,spineind,hbcomps(exptind));
                allheart=zeros(size(exptorder,1),length(megind),beatlen);
            else
                [heartest]=grb_est_heartbeat(D,spineind,hbcomps(exptind),beatlen,megind);
            end
            allheart(exptind,megind,:)=heartest;
        end % if HB


        %% downsample OPM data to match EMG

        S=[];
        S.fsample_new=EMGstruct.fsample;
        S.D=D;
        dD=spm_eeg_downsample(S);

        %% splice EMG data into the downsampled OPM dataset

        EMGchannames=strvcat(EMGchanname,otherEMGchan);
        [dD,EMGsamples,EMGdata]=addinEMG(D,dD,fullfile(EMGpath, EMGfilename),EMGtriglatency,OPMEMGsync,OPMchansreplace,EMGchannames);


        if EMGHP
            S=[];
            S.D=dD;
            S.band='high';
            S.order=5;
            S.freq=emghpf;
            S.chans2filter=dD.indchantype('EMG'); %
            dD=spm_eeg_filter(S);
        end 

        %% now make standard SPM datasets split into epochs of length EpochTime

        if ONEEPOCH %% SINGLE EPOCH PER DATASET
            epochsamples=length(EMGsamples);
        else
            epochsamples=EpochTime*dD.fsample;
        end

        Nepochs=floor(length(EMGsamples)/epochsamples);

        %make trl matrix for epoching
        trl=zeros(Nepochs,3);
        for f=1:Nepochs
            trl(f,:)=[(f-1)*epochsamples+min(EMGsamples) f*epochsamples+min(EMGsamples)-1 0];
        end

        S=[];
        S.D=dD;
        S.conditionlabels=[exptorder{exptind,2} task num2str(exptind)];
        S.trl=trl;
        S.bc=0; %
        S.prefix=sprintf('e%dms',round(EpochTime*1000));
        dDep = spm_eeg_epochs(S);




        %% remove outliers
        S=[];
        S.D=dDep;
        [dDep,retain] = spm_opm_removeOutlierTrials(S);
        savename=sprintf('retainedtrials%srun%s_%s',sub,run,pstr(1));
        save(fullfile(savedir,savename),'retain')


        allfilenames=strvcat(allfilenames,dDep.fname);

    end % for expt ind

    % now merge (before heartbeat removal)

    cd(dD.path);
    S=[];
    S.D=allfilenames;
    S.prefix=pstr;
    DAll=spm_eeg_merge(S);

         datft=spm2fieldtrip(DAll);
        datft=rmfield(datft,'hdr');
        cfg=[];
        cfg.channel=datft.grad.label;
        %dbfig=ft_databrowser(cfg,datft)

    if HB==1
        avheart=squeeze(mean(allheart));
        figure
        plot(avheart')

        BALANCE=1;
        DAll=grb_remove_heartbeat(DAll,avheart,spineind,megind,BALANCE);
    end

        datft=spm2fieldtrip(DAll);
        datft=rmfield(datft,'hdr');
        cfg=[];
        cfg.channel=datft.grad.label;
       % dbfig=ft_databrowser(cfg,datft)

    if HFCflag
        S=[];
        S.prefix='hfc';
        S.D=DAll;
        fprintf('\n Running HFC');
        [DAllhf, Yinds]=spm_opm_hfc(S);
        chans_hf = chanlabels(DAllhf,Yinds);

        S = [];
        S.D1=DAllhf;
        S.D2=DAll;
        S.plot = 1;
        S.channels = chans_hf;
        S.triallength = 2000;
        S.wind = @hanning;
        [shield,f] = spm_opm_rpsd(S);
        xlim([1,100])
        title('hfc shielding factor')
        DAll=DAllhf;

    else
        warning('no HFC')
    end
end % for RIGHT



%final merged file

