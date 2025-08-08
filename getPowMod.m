%% sensor level power modulations - brain and EMG
clear all
close all
clc

%subs = {'OP00212', 'OP00213',  'OP00215', 'OP00219', ...
         %'OP00221', 'OP00225', 'OP00226'};

subs={'OP00224'}; %_002

HFC=1;

for s=1:length(subs)
sub = subs{s};
save_dir = fullfile('D:\MSST001', [sub '_contrast']);

if HFC
filename = fullfile('D:\MSST001', ...
    ['sub-' sub], ...
    'ses-001', ...
    'meg', ...
    ['pmergedoe1000mspddfflo45hi45hfcstatic_002_array1.mat']);
else
warning('havent yet coded no hfc option')
end
DAll=spm_eeg_load(filename);


%%
load('D:\brainspineconnectivity\topoplot\brain_layout.mat')

%Zbrainchans=brainlay.label;

sensors=DAll.sensors('MEG');
alllabels=chanlabels(DAll);
headind=[];
for f = 1:length(sensors.chanpos)
chanind=find(strcmp(sensors.label(f),alllabels));
    if sensors.chanpos(f,2) > 200
        headind = [headind chanind];
    end
end

emgidx=find(contains(DAll.chanlabels,'EXG1'));

ftdat=spm2fieldtrip(DAll);

cfg = [];
cfg.channel = alllabels(headind(startsWith(alllabels(headind), 'Z')));
opmdat = ft_selectdata(cfg, ftdat);

cfg=[];
cfg.channel='EXG1';
cfg.rectify ='yes';
cfg.detrend='yes';
EMG=ft_preprocessing(cfg, ftdat);

alldat=ft_appenddata([], opmdat, EMG); 


% cfg            = [];
% cfg.output     = 'fourier';
% cfg.method     = 'mtmfft';
% cfg.foilim     = [10 35];
% cfg.tapsmofrq  = '2';
% cfg.keeptrials = 'yes';
% % cfg.pad='maxperlen';
% 
% freq_dat   = ft_freqanalysis(cfg, alldat);
% 
% cfg=[];
% stattrials=find(freq_dat.trialinfo==1);
% cfg.trials=stattrials;
% statdat=ft_selectdata(cfg,freq_dat);
% 
% seed='EXG1';
% brain=Zbrainchans;
% tmpcmb=repmat({seed}, length(brain),1);
% combs=[brain tmpcmb];
% 
% cfg            = [];
% cfg.method     = 'powcorr';
% cfg.channelcmb = combs;
% cfg.removemean='yes';
% coh      = ft_connectivityanalysis(cfg, statdat);
% 
% %figure; plot(coh.freq, coh.powcorrspctrm)
% 
% betafreqs=find(coh.freq >=15 & coh.freq<=30);
% betacoh=mean(coh.powcorrspctrm(:,betafreqs),2);
% 
% % 
% % figure
% % subplot(211)
% % topobrain(brainlay, betacoh)
% 
% % Sum coherence across frequencies for each channel
% coh_sum = sum(coh.powcorrspctrm, 2);  % sum over freqs (dim 2)
% 
% % Find (local) index of maximum sum
% [~, max_idx] = max(abs(coh_sum));
% brainlay.label(max_idx)
% 
% subplot(212)
% plot(coh.freq, coh.powcorrspctrm(max_idx,:))
% 
% 
% %% are there stat sig power modulations for brain chans and emg?
% 
% cfg            = [];
% cfg.output     = 'pow';
% cfg.method     = 'mtmfft';
% cfg.foilim     = [5 45];
% cfg.tapsmofrq  = '2';
% cfg.keeptrials = 'yes';
% 
% freq_dat2   = ft_freqanalysis(cfg, alldat);
% 
% %split data into two conditions 
% cfg=[];
% stattrials=find(freq_dat2.trialinfo==1);
% resttrials=find(freq_dat2.trialinfo==2);
% cfg.trials=stattrials;
% statdat=ft_selectdata(cfg,freq_dat2);
% cfg.trials=resttrials;
% restdat=ft_selectdata(cfg,freq_dat2);
% 
% 
% % average over freqs
% statpow=mean(statdat.powspctrm(:,:,:),3);
% restpow=mean(restdat.powspctrm(:,:,:),3);
% 
% ntrials=min([size(statpow,1) size(restpow,1)]);
% 
% statpow=log(statpow(1:ntrials,:));
% restpow=log(restpow(1:ntrials,:));
% 
% [~, p, ~, stats] = ttest2(statpow, restpow);
% warning('two-tailed t test')
% braint=stats.tstat(1:end-1);



% figure
% subplot(121)
% topobrain(brainlay, braint'); hold on
% for k=1:length(p)-1
% if p(k) < 0.05
% scatter(brainlay.pos(k,1), brainlay.pos(k,2), 'ro')
% end
% end
% title('Brain t-stats')
% 
% emg_stat = statpow(:, end); % Last channel is EMG
% emg_rest = restpow(:, end);
% emg_mean = [mean(emg_stat), mean(emg_rest)];
% emg_std = [std(emg_stat), std(emg_rest)] ./ sqrt(ntrials); % SEM

% Plot EMG power with p-value
% subplot(122)
% bar(emg_mean)
% hold on
% errorbar(1:2, emg_mean, emg_std, 'k.', 'LineWidth', 1.5)
% xticks([1 2])
% xticklabels({'Contraction', 'Rest'})
% ylabel('log Beta power')
% title(sprintf('EMG beta power (p = %.2f)', p(end))) % p-value for EMG
% 
% log_power_diff = mean(statpow, 1) - mean(restpow, 1);  % 1 x Nchannels
% 
% save(fullfile(save_dir, ['brainEMGmod_subject_' sub '.mat']), 'log_power_diff', 'betacoh');




%% CVA brain emg - spatial brai mixture

% filter to beta ish
cfg=[];
cfg.bpfilter='yes';
cfg.bpfreq=[10 35];
EMG=ft_preprocessing(cfg,EMG);
opmdat=ft_preprocessing(cfg,opmdat);

EMGdat=[];
opmdat1=[];
for t=1:length(EMG.trial)
    EMGdat(:,t)=EMG.trial{t}(:,:);
    opmdat1(:,:,t)=opmdat.trial{t}(:,:);
end

Y=squeeze(mean(opmdat1,2));
X=squeeze(mean(EMGdat,1));

Y=(Y-mean(Y))';
X=(X-mean(X))';
CVA=spm_cva(Y,X); 

fprintf('CVA p is %g for subject %s\n', CVA.p,sub)

end
