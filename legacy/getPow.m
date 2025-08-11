function [trialPow, chanlabs2use] = getPow(DAll,labels)

%% sensor level coherence
%brainchans=labelToIndex(labels);
brainchanidx=find(contains(DAll.chanlabels,labels));
brainchanidx = setdiff(brainchanidx, DAll.badchannels);

brain_labels=DAll.chanlabels(brainchanidx);
emgidx=find(contains(DAll.chanlabels,'EXG1'));


ftdat=spm2fieldtrip(DAll);

% could add excluding bad channels here

cfg=[];
cfg.channel=ftdat.label(~contains(ftdat.label,'EXG1'));
opmdat=ft_selectdata(cfg, ftdat);

cfg=[];
cfg.channel='EXG1';
cfg.rectify ='yes';
cfg.detrend='yes';
EMG=ft_preprocessing(cfg, ftdat);

alldat=ft_appenddata([], opmdat, EMG); %use this for Neurospec
chanlabs2use=brain_labels(contains(brain_labels, 'Z')); %radial
    
    cfg            = [];
    cfg.output     = 'pow';
    cfg.channel= chanlabs2use;
    cfg.method     = 'mtmfft';
    cfg.foilim     = [10 35];
    cfg.tapsmofrq  = '2';
    cfg.keeptrials = 'yes';
   % cfg.pad='maxperlen';

    freq_dat   = ft_freqanalysis(cfg, alldat);

    trialPow=squeeze(mean(freq_dat.powspctrm,3)); %trials by chans

%     seed='EXG1';
%     brain=alldat.label(contains(alldat.label,labels));
%     tmpcmb=repmat({seed}, length(brain),1);
%     combs=[brain tmpcmb];
% 
%     cfg            = [];
%     cfg.method     = 'coh';
%     cfg.complex    ='absimag';
%     cfg.channelcmb = combs;
%     coh      = ft_connectivityanalysis(cfg, freq_dat);
    
%     fig=figure;
%     cfg=[];
%     ft_connectivityplot(cfg,coh)
%     waitfor(fig)

%     for k=1:length(brain)
%     f=figure;
%     plot(coh.freq, coh.cohspctrm(k,:))
%     waitfor(f)
%     end


% figure; plot(coh.freq, coh.cohspctrm)
% 
% % Sum coherence across frequencies for each channel
% coh_sum = sum(coh.cohspctrm, 2);  % sum over freqs (dim 2)
% 
% % Find (local) index of maximum sum
% [~, max_idx] = max(coh_sum);
% 
% brainchan2use=labels(max_idx);
% fprintf('brainchan to use:%s/n', brainchan2use)

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