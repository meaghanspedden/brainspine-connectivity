%% process data from hand contraction. ~ July 2023
% go through either all right or all left contractions preprocess and save
%R or L--> line 27

clear all;
close all;

sub='OP00116';
bids_session='001';
dsprefix=0;
brainchan_labels={'19','1B', 'MW', 'MY','OK','MZ','35','DI','MU','DQ','OH','DG','A1'};
badchans={};
posfileneck='D:\MSST001\sub-OP00116\ses-001\meg\sub-OP00116_ses-001_positions_new.tsv';
EMGpath='D:\OP00116_experiment\EMGfiles';
EMGfiletemplate='000116_static_';
savedir='D:\MSST001\Coh_results00116';

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

%% which opm runs correspond to which emg runs

Bothexptorder={...
    '001' 'R' '07' ;
    '002'  'R'  '09';
    '003'  'R'  '11';
    '004'  'R'  '13';
    '001'  'L'  '23';
    '002'  'L'  '25';
    '003'  'L'  '27'};  %% 116 left only has 3 runs

%% trig info

OPMEMGsync='NI-TRIG-1'; %% EMG trigger syncing to OPM
OPMchansreplace=strvcat('NI-TRIG-2','NI-TRIG-3'); %% replace this/these OPM channels with EMG data
EMGtriglatency=0.05;

%% paths
datadir= 'D:\MSST001';
addpath('D:\spm')
spm('defaults','EEG')
addpath(genpath('D:\brainspineconnectivity'))
%warning('MAKE SURE TO CHANGE spm_eeg_filter to GRB version')


%% right and left hands

for RIGHT=1

    if RIGHT
        exptorder=Bothexptorder(1:4,:);
        pstr='rhand';
    else
        exptorder=Bothexptorder(5:7,:);
        pstr='lhand';
    end

    %%

    allfilenames=[];allxfilenames=[];
    for exptind=1:size(exptorder,1) %loop through runs

        if strcmp(exptorder(exptind,2),'R')
            EMGchanname='EMG 1';
            otherEMGchan='EMG 2';
            EMGfilename=[EMGfiletemplate,cell2mat(exptorder(exptind,3)),'.smrx'];
            task='staticright';
            hbcomp=[1 1 1 1];
        end

        if strcmp(exptorder(exptind,2),'L')
            EMGchanname='EMG 2';
            otherEMGchan='EMG 1';
            EMGfilename=[EMGfiletemplate,cell2mat(exptorder(exptind,3)),'.smrx'];
            task='staticleft';
            hbcomp=[1 1 1];
        end

        run=cell2mat(exptorder(exptind,1));

        %% make SPM object
        [D,fullfilename]=convert_opm2spm(datadir,posfileneck,sub,bids_session,task,run,'single',dsprefix);

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
        title('raw psd')

        
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

        spineind=D.indchannel(D.sensors('MEG').label);
        allspineind(exptind,:)=spineind;
        allchannames(exptind,:)=D.sensors('MEG').label;


        %% mark bad channels

        if~isempty(badchans)
            badind=[];

            for f=1:size(badchans,1)
                badind(f)=find(contains(D.chanlabels,deblank(badchans(f,:))));
            end
            D = badchannels(D, badind, 1); %% set channels to bad
        end

        %%note that I have modified this in grb_est_heartbeat

        if HB % estimate heartbeat over all channels (will remove after merging files)
            if exptind==1
                [heartest,beatlen,megind]=grb_est_heartbeat(D,spineind,hbcomp(exptind));
                allheart=zeros(size(exptorder,1),length(megind),beatlen);
            else
                [heartest]=grb_est_heartbeat(D,spineind,hbcomp(exptind),beatlen,megind);
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
        end; % if EMG HP



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
        [dDep, retain] = spm_opm_removeOutlierTrials(S);
        
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

    if HB==1
        avheart=squeeze(mean(allheart));
        figure
        plot(avheart')

        BALANCE=1;
        Dbeforeproject=DAll;
        DAll=grb_remove_heartbeat(DAll,avheart,spineind,megind,BALANCE);
    end



    if HFCflag
        S=[];
        S.prefix='hfc';
        S.D=DAll;
        fprintf('\n Running HFC');
        [DAllhf, Yinds]=spm_opm_hfc(S);
        chans_hf = chanlabels(DAllhf,Yinds);

        S = [];
        S.D1 = DAllhf;
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

