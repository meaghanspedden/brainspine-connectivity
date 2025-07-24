%% sensor level coherence


sub = 'OP00215';
concat = 0;  % Set to 0 to load only one run
run2use='001';

if concat %note that this has hb removed and has rest!!!

    filename = fullfile('D:\MSST001', ...
        ['sub-' sub], ...
        'ses-001', ...
        'meg', ...
        ['pmergedoe1000mspddfflo45hi45hfcstatic_001_array1.mat']);

    DAll=spm_eeg_load(filename);

else

    filename = fullfile('D:\MSST001', ...
        ['sub-' sub], ...
        'ses-001', ...
        'meg', ...
        ['oe1000mspddfflo45hi45hfcstatic_' run2use '_array1.mat']);

    DAll=spm_eeg_load(filename);
end

%%

brainchan_labels={'H1', 'H6', 'B3', 'H8', 'G3', 'B6', 'H7', 'H4', 'G4', 'B4', 'B2', 'H2', 'B7', 'B8', 'B5', 'B1'};

brainchans=labelToIndex(brainchan_labels);
brainchanidx=find(contains(DAll.chanlabels,brainchans));
emgidx=find(contains(DAll.chanlabels,'EXG1'));

ftdat=spm2fieldtrip(DAll);
if concat %only static condition
    cfg=[];
    cfg.trials=find(ftdat.trialinfo==1);
    ftdat=ft_selectdata(cfg,ftdat);
end

% could add excluding bad channels here
Zchans=ftdat.label(contains(ftdat.label,brainchans) & contains(ftdat.label, 'Z'));

cfg=[];
cfg.channel=Zchans;
opmdat=ft_selectdata(cfg, ftdat);

cfg=[];
cfg.channel='EXG1';
cfg.rectify ='yes';
cfg.detrend='yes';
EMG=ft_preprocessing(cfg, ftdat);

alldat=ft_appenddata([], opmdat, EMG); %use this for Neurospec


cfg            = [];
cfg.output     = 'fourier';
cfg.method     = 'mtmfft';
cfg.foilim     = [5 45];
cfg.tapsmofrq  = '2';
cfg.keeptrials = 'yes';
% cfg.pad='maxperlen';

freq_dat   = ft_freqanalysis(cfg, alldat);



seed='EXG1';
brain=Zchans;
tmpcmb=repmat({seed}, length(brain),1);
combs=[brain tmpcmb];

cfg            = [];
cfg.method     = 'coh';
cfg.complex    ='absimag';
cfg.channelcmb = combs;
coh      = ft_connectivityanalysis(cfg, freq_dat);

%     fig=figure;
%     cfg=[];
%     ft_connectivityplot(cfg,coh)
%     waitfor(fig)

% for k=1:length(brain)
%     f=figure;
%     plot(coh.freq, coh.cohspctrm(k,:),'LineWidth',2)
%     title(sprintf('%s',Zchans{k}))
%     waitfor(f)
% end


figure; plot(coh.freq, coh.cohspctrm)

% Sum coherence across frequencies for each channel
coh_sum = sum(coh.cohspctrm, 2);  % sum over freqs (dim 2)

% Find (local) index of maximum sum
[~, max_idx] = max(coh_sum);

% brainchan2use=labels(max_idx);
% fprintf('brainchan max coh:%s/n', brainchan2use)

%neurospec coh

% addpath('D:\neurospec20')
% addpath(genpath('D:\neurospec211NEW'))
%
% contdat=[];
% for k=1:length(alldat.trial)
%     contdat=[contdat alldat.trial{k}];
% end
%
%
% opt_str='t2 M5';
% seg_pwr=10;
%
% cohmat=[];
% f11=[];
% f22=[];
%
% EMG=contdat(emgidx,:);
%
% for cc=1:length(brainchanidx)
%     dat1=contdat(brainchanidx(cc),:)';
%     [f,t,cl,sc] = sp2a2_R2_mt(dat1,EMG',alldat.fsample,seg_pwr,opt_str);
%     cohmat(:,cc)=f(:,11)+f(:,12);
%
%     %f(:,4);
%     f11(:,cc) = f(:,2);
%     f22(:,cc) = f(:,3);
% end
%
% % for k=1:length(brainchanidx)
% %     fig=figure;
% %     plot(f(:,1), cohmat(:,k))
% %     xlim([0 100])
% %     waitfor(fig)
% % end
%
% medcoh=median(cohmat,2);
%
% figure; plot(f(:,1), cohmat)
% hold on
% plot(f(:,1), medcoh, 'k', 'LineWidth',2)
% xlim([0 100])