%% coherence

clc
%clear all
close all

D=spm_eeg_load('D:\MSST001\sub-OP00189\ses-001\meg\mergeoe1000mse1000msdfflo45hi45sub-OP00189_ses-001_task-static_run-001_meg.mat');

%brainchan_labels
%ref_labels
%sc_labels

emg_label=all_labels(EMGidx);

goodchanidx=setdiff(1:size(D,1),badchannels(D));

ftdat=spm2fieldtrip(D);
ftdat=rmfield(ftdat,'hdr');

cfg=[];
cfg.channel=D.chanlabels(goodchanidx);
ft_selectdata(cfg,ftdat);

%% freq

cfg            = [];
cfg.output     = 'fourier';
cfg.method     = 'mtmfft';
cfg.foilim     = [0 200];
cfg.tapsmofrq  = 2;
cfg.keeptrials = 'yes';

freq_dat   = ft_freqanalysis(cfg, ftdat);

%% BRAIN-----------------------------------------------

seed=emg_label;
megchans=brainchan_labels';
tmpcmb=repmat(seed, length(megchans),1);
combs=[megchans tmpcmb];

cfg            = [];
cfg.method     = 'coh';
cfg.complex    ='absimag';
cfg.channelcmb = combs;
coh      = ft_connectivityanalysis(cfg, freq_dat);

%%

lay = make_layout_brainchans(ftdat,brainchan_labels);


fig=figure;
fig.Position = [100, 200, 600, 400];
cfg                  = [];
cfg.parameter        = 'cohspctrm';
cfg.xlim             = [15 25];
cfg.zlim             ='zeromax';
cfg.refchannel       = seed;
cfg.interplimits     ='electrodes';
cfg.layout           = lay;

ft_topoplotER(cfg,coh)
colorbar


maxchan='11-B3-Z';
maxchanidx=find(contains(coh.labelcmb,maxchan));

figure;
plot(coh.freq, coh.cohspctrm(maxchanidx,:), 'LineWidth', 2, 'Color', [0, 0.5, 0.4]);
box off;
xlim([0 60]);
xlabel('Frequency (Hz)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Coherence', 'FontSize', 14, 'FontWeight', 'bold');
title('Coherence Spectrum', 'FontSize', 16, 'FontWeight', 'bold');
grid on;
set(gca, 'GridColor', [0.7, 0.7, 0.7]);
set(gca, 'GridLineStyle', '--');
set(gca, 'FontSize', 12, 'FontWeight', 'bold');
xticks(0:10:60);
yticks(0:0.1:1);


%% SPINAL CORD

seed=emg_label;
megchans=sc_labels';
tmpcmb=repmat(seed, length(megchans),1);
combs=[megchans tmpcmb];

cfg            = [];
cfg.method     = 'coh';
cfg.complex    ='absimag';
cfg.channelcmb = combs;
coh      = ft_connectivityanalysis(cfg, freq_dat);

%%

lay_sc = make_layout_neckchans(ftdat,sc_labels, '-Y');


fig=figure;
fig.Position = [100, 200, 600, 400];
cfg                  = [];
cfg.parameter        = 'cohspctrm';
cfg.xlim             = [5 45];
cfg.zlim             ='zeromax';
cfg.refchannel       = seed;
cfg.interplimits     ='electrodes';
cfg.layout           = lay_sc;
cfg.marker = 'labels';
% 
% ft_topoplotER(cfg,coh)
% colorbar

cfg.layout='ordered';
ft_multiplotER(cfg,coh)



% maxchan='17-C1-Y';
% maxchanidx=find(contains(coh.labelcmb,maxchan));
% 
% figure;
% plot(coh.freq, coh.cohspctrm(maxchanidx,:), 'LineWidth', 2, 'Color', [0, 0.5, 0.4]);
% box off;
% xlim([0 60]);
% xlabel('Frequency (Hz)', 'FontSize', 14, 'FontWeight', 'bold');
% ylabel('Coherence', 'FontSize', 14, 'FontWeight', 'bold');
% title('Coherence Spectrum', 'FontSize', 16, 'FontWeight', 'bold');
% grid on;
% set(gca, 'GridColor', [0.7, 0.7, 0.7]);
% set(gca, 'GridLineStyle', '--');
% set(gca, 'FontSize', 12, 'FontWeight', 'bold');
% xticks(0:10:60);
% yticks(0:0.1:1);


%%

maxchan='-Z';
maxchanidx=find(contains(coh.labelcmb,maxchan));

figure;
plot(coh.freq, coh.cohspctrm(maxchanidx,:), 'LineWidth', 2, 'Color', [0, 0.5, 0.4]);
box off;
xlim([0 40]);
xlabel('Frequency (Hz)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Coherence', 'FontSize', 14, 'FontWeight', 'bold');
title('Coherence Spectrum', 'FontSize', 16, 'FontWeight', 'bold');
grid on;
set(gca, 'GridColor', [0.7, 0.7, 0.7]);
set(gca, 'GridLineStyle', '--');
set(gca, 'FontSize', 12, 'FontWeight', 'bold');
xticks(0:10:60);
%yticks(0:0.1:1);
