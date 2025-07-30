%% opm static contraction analysis

clear all
close all

addpath('D:\spm')
spm('defaults','EEG')

addpath(genpath('D:\brainspineconnectivity\plotting'))

%MEGruns={'001'};

refsens={'X51', 'Y51', 'Z51','X49', 'Y49','Z49','X53', 'Y53', 'Z53'};

headsens_col={'R8', 'R4', 'R1', 'R2', 'R3','R5', 'R6', 'R7'};

spinesens_col={'DB2', 'DB4', 'DB6', 'DB7', 'DB8','BR1', 'BR2', 'BR3',...
    'BR4', 'BR5', 'BR6', 'BR7', 'BR8', 'LB1', 'LB2', 'LB3', 'LB4',...
    'LB5', 'LB6', 'LB7', 'LB8', 'G1', 'G2' 'G3', 'G4', 'G5', 'G6', 'G7', 'G8',...
    'SW1', 'SW2', 'SW3', 'SW4', 'SW5', 'SW6', 'SW7', 'SW8'};

badchans={'X9', 'Y10', 'Y48', 'X16', 'Y41', 'Y16', 'X48','Z41', 'X10',...
    'Z48', 'Z16', 'Z10', 'X41','Y28','X28', 'Z28','X11', 'Y56', 'Z13', 'Z12',...
    'Y22', 'X22', 'Y13', 'X14', 'Y14', 'X12','Y12', 'Y9', 'Z15', 'X15', 'Y15',...
    'Y18', 'Y30', 'X13', 'X18', 'Z14'};

headsens={};
for ii=1:length(headsens_col)
    sensornumber = get_sensornumber(headsens_col{ii});
    headsens=[headsens sensornumber];
end


spinesens={};
for ii=1:length(spinesens_col)
    sensornumber = get_sensornumber(spinesens_col{ii});
    spinesens=[spinesens sensornumber];
end
%%
%will need to add positions here
S = [];
S.data = 'Static-run-001_array1.lvm';
D = spm_opm_create(S);

labs=chanlabels(D);
badchanidx=find(contains(labs,badchans));

EMGchan='A8';
EMGidx=find(strcmp(labs, EMGchan));

%figure; plot(D.time, D(EMGidx,:))

%%

MEGchans=find(contains(D.chantype,'MEGMAG')); %idx
chanidxtoplot =setdiff(MEGchans,badchanidx); %includes ref sens

%%
% look at psd
S = [];
S.D = D;
S.channels = D.chanlabels(chanidxtoplot);
S.plot = 1;
S.triallength = 2000;
S.wind = @hanning;
spm_opm_psd(S);
xlim([1,100])
title('raw')


%% high pass
S = [];
S.D = D;
S.type = 'butterworth';
S.band = 'high';
S.freq = 5;
S.dir = 'twopass';
Dfilt = spm_eeg_filter(S);
%% low pass
S = [];
S.D = Dfilt;
S.type = 'butterworth';
S.band = 'low';
S.freq = 45;
S.dir = 'twopass';
Dfilt = spm_eeg_filter(S);

%% band stop

S = [];
S.D = Dfilt;
S.type = 'butterworth';
S.band = 'stop';
S.freq = [49 51];
S.dir = 'twopass';
Dfilt = spm_eeg_filter(S);

%% need to use this to further identify bad channels...

% ftdat=spm2fieldtrip(Dfilt);
% ftdat=rmfield(ftdat,'hdr');
% 
% cfg=[];
% cfg.channel=ftdat.label(chanidxtoplot);
% ft_databrowser(cfg,ftdat)

%% synthetic gradiometry

refidx=find(contains(Dfilt.chanlabels,refsens));
%Dfilt=chantype(Dfilt,refidx,'REF'); %if you want to use ref sens not 1 pc

ref_ts=Dfilt(refidx,:); %time series for 9 ref channels
ref_ts=ref_ts-mean(ref_ts,2); %mean centre
ref_cv=ref_ts*ref_ts'; %cov matrix c x c

[U,~,~]=svd(ref_cv);

replacedat=U(:,1); %% first pC
replaceind=find(contains(labs,'Data_Valid2')); %channel to replace

ftdat=spm2fieldtrip(Dfilt);
ftdat.trial{1}(replaceind,:)=replacedat'*ref_ts;
ftdat.label(replaceind)={'refPC1'};
ftdat=rmfield(ftdat, "hdr");

Dfilt=spm_eeg_ft2spm(ftdat, 'testconv');
refidx1=find(contains(Dfilt.chanlabels,'refPC1'));
Dfilt=chantype(Dfilt,refidx1,'REF');
Dfilt=chantype(Dfilt, MEGchans, 'MEGMAG');
Dfilt=chantype(Dfilt, refidx, 'Other');

S = [];
S.D = Dfilt;
S.plot = 1;
S.channels = Dfilt.chanlabels(chanidxtoplot);
S.triallength = 2000;
S.wind = @hanning;
spm_opm_psd(S);
xlim([1,100])
title('pre SG')



S=[];
S.D=Dfilt;
S.confounds={'REF'};
D_SG = spm_opm_synth_gradiometer(S);


S = [];
S.D = D_SG;
S.plot = 1;
S.channels = D_SG.chanlabels(chanidxtoplot);
S.triallength = 2000;
S.wind = @hanning;
spm_opm_psd(S);
xlim([1,100])
title('post SG')

S = [];
S.D1 = Dfilt;
S.D2=D_SG;
S.plot = 1;
S.channels = D_SG.chanlabels(chanidxtoplot);
S.triallength = 2000;
S.wind = @hanning;
[shield,f] = spm_opm_rpsd(S);

xlim([1,100])


%% hb removal

[heartep,beatlen]=mes_est_heartbeat(ftdat,hbcomp,beatlen)

