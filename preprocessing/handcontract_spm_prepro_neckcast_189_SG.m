% new coherence pipeline Feb 2025

%need a good way to sort how to deal with bad positions? could update pos
%file with 
clear all;
close all;

sub='OP00189';
bids_session='001';
task='static';

brainchan_labels={'13-B5', '15-B7', '9-B1', '11-B3', '14-B6',...
    '32-D8', '28-D4', '25-D1', '27-D3', '31-D7', '29-D5', '30-D6', '26-D2'};
badchans={'B8', 'D4', 'F8', 'B2', 'F1', 'B1', 'B4', 'G3', '58-H2-Z','60-H4-Z','61-H5-Z', '62-H6-Z', '63-H7-Z', '63-H7-Z',...
    '64-H8-Z'};


ref_labels={'G5', 'G1'};
EMGchan='A8';

hbcomp=[1 1 1 1];

posfile='D:\MSST001\sub-OP00189\ses-001\meg\positions.tsv';
savedir='D:\MSST001\Coh_results00189';

if ~exist(savedir, 'dir')
    mkdir(savedir);
end

%% Processing Options

fband=[5 45]; %% OPM filter band
stpband1=[48 52];
stpband2=[98 102];
EpochTime=1; %  %% length of epoch in seconds

% flags to turn these filters on/off
HFCflag=0; %% homogenous field correction (0 or 1)
MEGHP=1;
MEGLP=1;
MEGBS=1; % stop band
HB=0; %% heartbeat removal
SG=1;

%% paths
datadir= 'D:\MSST001';
addpath('D:\spm')
spm('defaults','EEG')
addpath(genpath('D:\brainspineconnectivity'))

exptruns={'001', '002', '003', '004', '005'};

allfilenames=[];


for thisrun=1:length(exptruns)

    close all
    run=exptruns{thisrun};

    %% make SPM object
    [D,fullfilename]=convert_opm2spm(datadir,posfile,sub,bids_session,task,run,'single');


    % get labels and indices for ref, spinal cord, and brain channels
    all_labels=D.chanlabels;
    chan_types=chantype(D);
    megmag_indices = strcmp(chan_types, 'MEGMAG'); %all opms including ref.

    opm_labels=all_labels(megmag_indices);
    brain_idx=find(contains(all_labels, brainchan_labels));
    sc_labels=opm_labels(~contains(opm_labels, brainchan_labels) & ~contains(opm_labels,ref_labels));
    brainchan_labels=all_labels(brain_idx);

    sc_idx=find(contains(all_labels,sc_labels));
    brain_idx=find(contains(all_labels, brainchan_labels));
    ref_idx=find(contains(D.chanlabels,ref_labels));

    %% psd

    S = [];
    S.D = D;
    S.plot = 1;
    S.channels = opm_labels;
    S.triallength = 2000;
    S.wind = @hanning;
    [po,freq] = spm_opm_psd(S);
    xlim([1,100])
    title('raw psd')


    % Find bad channels based on psd
    freqTarget = 40;
    freqIdx = dsearchn(freq(:), freqTarget); 
    powerThresholdHigh = 10^2.5;
    powerThresholdLow = 10^1;

    badchanIndices = find(po(freqIdx, :) > powerThresholdHigh | po(freqIdx, :) < powerThresholdLow);
    
    %global indices
    badchanIdx1 = find(contains(all_labels, badchanLabels)); 
  
   % Identify bad channels noted during the experiment
    badchanIdx2 = find(contains(all_labels, badchans));

    allBadchanIndices = unique([badchanIndices1, badchanIndices2]);
    badchanLabels = opm_labels(allBadchanIndices); 
    D = badchannels(D, allBadchanIndices, 1);

    % Determine the good channels
    goodopmLabels = setdiff(opm_labels, badchanLabels);
    goodopmIdx = find(contains(all_labels, goodopmLabels));  

    %update good spinal cord and brain channel labels
    sc_labels=goodopmLabels(contains(goodopmLabels,sc_labels));
    brainchan_labels=goodopmLabels(contains(goodopmLabels, brainchan_labels));
    

    S = [];
    S.D = D;
    S.plot = 1;
    S.channels = goodopmLabels;
    S.triallength = 2000;
    S.wind = @hanning;
    [po,freq] = spm_opm_psd(S);
    xlim([1,100])
    title('raw psd')


    %find EMG channel and change type to EMG
    %
    % for i = 1:length(D.chanlabels)
    %     if contains(D.chanlabels{i}, 'A')  % Check if 'A' is in the channel label
    %         f=figure;  % Create a new figure for each plot
    %         plot(D.time, D(i,:));  % Plot the data for the corresponding channel
    %         title(['Channel: ', D.chanlabels{i}]);  % Title with channel label
    %         xlabel('Time (s)');
    %         ylabel('Signal Value');
    %         waitfor(f)
    %     end
    % end

    %channel for EMG is A8

    EMGidx=find(contains(D.chanlabels, EMGchan));
    D = chantype(D,EMGidx, 'EMG'); %have checked that spm filters channels labelled as EMG

    if MEGHP
        S=[];
        S.D=D;
        S.freq=min(fband);
        S.order=5;
        S.band='high';
        S.prefix=sprintf('hi%d',fband(2));
        D=spm_eeg_filter(S);
    end


    %% visual inspection

    ftdat=spm2fieldtrip(D);
    ftdat=rmfield(ftdat,'hdr');

    cfg=[];
    cfg.channel=goodopmLabels;
    ft_databrowser(cfg,ftdat)

    % identify any bad channels and add to bad chan list...




    %%  NOW LP FILTER
    if MEGLP
        S=[];
        S.D=D;
        S.freq=max(fband);
        S.order=5;
        S.band='low';
        S.prefix=sprintf('lo%d',fband(2));
        D=spm_eeg_filter(S);
    end


    %% BAND STOP
    if MEGBS
        S=[];
        S.D=D;
        S.band='stop';
        S.order=3;
        S.freq=stpband1;
        S.D=D;
        D=spm_eeg_filter(S);

        S=[];
        S.D=D;
        S.band='stop';
        S.order=3;
        S.freq=stpband2;
        D=spm_eeg_filter(S);

    end



    if HB % estimate heartbeat over all channels (will remove after merging files)
        if thisrun==1
            [heartest,beatlen,megind]=grb_est_heartbeat(D,sc_idx,hbcomp(thisrun), [],brain_idx);
            allheart=zeros(size(exptorder,1),length(megind),beatlen);
        else
            [heartest]=grb_est_heartbeat(D,spineind,hbcomp(exptind),beatlen,megind);
        end
        allheart(exptind,megind,:)=heartest;
    end % if HB



    %% synthetic gradiometry

    if SG==1

        D=chantype(D,refidx,'REF');

        S=[];
        S.D=D;
        S.confounds={'REF'};

        D_SG = spm_opm_synth_gradiometer(S);

        S = [];
        S.D1 = D;
        S.D2=D_SG;
        S.plot = 1;
        S.channels = goodopmLabels;
        S.triallength = 2000;
        S.wind = @hanning;
        [shield,f] = spm_opm_rpsd(S);
        xlim([1,100])
        title('shielding factor (dB)')

        D=D_SG; 
    end


    %% now make SPM datasets split into epochs of length 1 s

    epoch_length = EpochTime * D.fsample; %samples
    N = length(D.time); % Total number of samples

    % Generate trial start and end indices
    trl = [(1:epoch_length:N-epoch_length)' (epoch_length:epoch_length:N)' zeros(floor(N/epoch_length),1)];

    S=[];
    S.D=D;
    S.trl=trl;
    S.bc=0; 
    S.prefix=sprintf('e%dms',round(EpochTime*1000));
    Depoch = spm_eeg_epochs(S);


    %% remove first 4 trials (EMG time to stabilise) and remove outliers

    [a,b,c]= fileparts(fullfile(Depoch));
    fname = fullfile(a,[S.prefix,b,c]);

    S=[];
    S.D=Depoch;
    S.inds = 5:size(Depoch,3);
    S.fname= fname;
    Depoch = spm_opm_selectTrials(S);

    %now remove outliers (opm data)
    S=[];
    S.D=Depoch;
    [Depoch, retain] = spm_opm_removeOutlierTrials(S);

    allfilenames=strvcat(allfilenames,Depoch.fname);

    ftdat=spm2fieldtrip(Depoch);
    ftdat=rmfield(ftdat,'hdr');

    cfg=[];
    cfg.channel=goodopmLabels;
    ft_databrowser(cfg,ftdat)

end % for opm runs

% now merge runs

cd(Depoch.path);

S=[];
S.D=allfilenames;
S.prefix='merge';
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
    S.D1 =DAllhf;
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


%final visual inspection--------------

ftdat=spm2fieldtrip(DAll);
ftdat=rmfield(ftdat,'hdr');

cfg=[];
cfg.channel=DAll.chanlabels(EMGidx);
ft_databrowser(cfg,ftdat)


%rectify EMG
DAll(EMGidx,:,:) =abs(DAll(EMGidx,:,:));
save(DAll)

