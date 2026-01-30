%% load all three time series and look at coherence between them
clear all
close all

subs = {'OP00212','OP00213','OP00215', 'OP00219', ...
 'OP00225', 'OP00221', 'OP00224'};

save_dir='C:\Users\mspedden\Documents\brainspine_save\';

addpath(genpath('C:\Users\mspedden\Documents\neurospec211NEW\neurospec211'))
addpath('C:\Users\mspedden\Documents\spm')
spm('defaults','EEG')

addpath('C:\Users\mspedden\Documents\fieldtrip')
ft_defaults

results=struct();

for ss=1:length(subs)

    sub=subs{ss};

    if strcmp(sub, 'OP00224')
        datwithEMGmerged = fullfile('C:\Users\mspedden\Documents', ...
            ['sub-' sub], ...
            'ses-001', ...
            'meg', ...
            'pmergedoe1000mspddfflo45hi45hfcstatic_002_array1.mat');
    else
        datwithEMGmerged = fullfile('C:\Users\mspedden\Documents', ...
            ['sub-' sub], ...
            'ses-001', ...
            'meg', ...
            'pmergedoe1000mspddfflo45hi45hfcstatic_001_array1.mat');
    end

 D=spm_eeg_load(datwithEMGmerged);


%load VE signals
filesp=sprintf('VE_spine_sub%s_forspectraVE_spine_subVE_forspectra.mat', sub);
sVE=load(fullfile(save_dir,filesp));
filebr=sprintf('sub%s_VE_brain_forspectra.mat', sub);
bVE=load(fullfile(save_dir,filebr));

bVE=bVE.VE;
sVE=sVE.VE;

ftdat = spm2fieldtrip(D);

%select EMG
cfg=[];
cfg.channel='EXG1';
EMG=ft_selectdata(cfg,ftdat);


cfg=[];
cfg.rectify='yes';
EMG=ft_preprocessing(cfg,EMG);


%% separate conditions

 statidx=find(ftdat.trialinfo==1);
 restidx=find(ftdat.trialinfo==2);
 [nTrials,~] = min([length(statidx) length(restidx)]);

 cfg=[];
 cfg.trials=statidx(1:nTrials);
 statB=ft_selectdata(cfg,bVE);
 statS=ft_selectdata(cfg,sVE);
 statEMG=ft_selectdata(cfg,EMG); 

 cfg.trials=restidx(1:nTrials);
 restEMG=ft_selectdata(cfg,EMG);
 restB=ft_selectdata(cfg,bVE);
 restS=ft_selectdata(cfg,sVE);

 %% fieldtrip functional connectivity
statB.label{1}='brain';
statS.label{1}='spine';
statEMG.label{1}='EMG';

alldat=ft_appenddata([],statB,statS,statEMG);

cfg            = [];
cfg.output     = 'fourier';
cfg.method     = 'mtmfft';
cfg.foilim     = [2 75];
cfg.tapsmofrq  = 2;
cfg.keeptrials = 'yes';

freq    = ft_freqanalysis(cfg, alldat);

cfg            = [];
cfg.method     = 'coh';
coh      = ft_connectivityanalysis(cfg, freq);


cfg=[];
cfg.method='coh';
cfg.labelcmb={'brain', 'EMG'};
cfg.partchannel={'spine'};
cohpart=ft_connectivityanalysis(cfg,freq);

figure; plot(coh.freq, squeeze(coh.cohspctrm(1,3,:)))
hold on
plot(cohpart.freq, squeeze(cohpart.cohspctrm(1,2,:)))
legend({'full', 'partial'})
xlim([10 35])


%% neurospec
 statBcont=[statB.trial{:}]; %format continuous
 statScont=[statS.trial{:}];
 statEMGcont=[statEMG.trial{:}];

% rectify
statEMGcont=abs(statEMGcont);

samp_rate=ftdat.hdr.Fs;
seg_pwr=11; 
opt_str='M3';
[f1,t1,cl1]=sp2a2_R2_mt(statBcont',statEMGcont',samp_rate,seg_pwr);
[f2,t2,cl2]=sp2a2_R2_mt(statBcont',statScont',samp_rate,seg_pwr,opt_str);
[f3,t3,cl3]=sp2a2_R2_mt(statScont',statEMGcont',samp_rate,seg_pwr,opt_str);


[fp tp clp] = sp2a2_R2_pc1(statBcont',statEMGcont',statScont',samp_rate,seg_pwr);

figure; plot(fp(:,1), fp(:,4)); hold on
plot(f1(:,1), f1(:,4))
legend({'partial', 'full'})
xlim([0 45])
% Example single-subject plot
freq_band = [10 35];  % Hz

figure('Color','w','Position',[100 100 600 800]);  % taller for vertical stacking

% Helper function to plot a line with color only in-band
plot_band = @(f, y, band, clr) ...
    plot(f, y, 'Color', [0.8 0.8 0.8], 'LineWidth', 2);  % plot all grey first

nSub = length(f1); % length of frequency vector

plot_colored_band = @(f, y, band, clr) ...
    plot(f(f>=band(1) & f<=band(2)), y(f>=band(1) & f<=band(2)), 'Color', clr, 'LineWidth', 2);

% -----------------------
% Brain ↔ EMG
subplot(3,1,1); hold on;

plot_band(f1(:,1), f1(:,12), freq_band, 'b');  % Reverse in grey first
plot_band(f1(:,1), f1(:,11), freq_band, 'r');  % Forward in grey

hRev = plot_colored_band(f1(:,1), f1(:,12), freq_band, 'b');  % Reverse in-band
hFwd = plot_colored_band(f1(:,1), f1(:,11), freq_band, 'r');  % Forward in-band

xlim([0 40]);
ylim([0 max([f1(:,11); f1(:,12)])*1.1]);

xlabel('Frequency (Hz)','FontSize',14);
ylabel('Coherence','FontSize',14);
title('Brain ↔ EMG','FontSize',16);
box off; set(gca,'FontSize',12);
legend([hRev hFwd],{'Reverse','Forward'},'Location','northeast');

% -----------------------
% Brain ↔ Spine
subplot(3,1,2); hold on;

plot_band(f2(:,1), f2(:,12), freq_band, 'b');
plot_band(f2(:,1), f2(:,11), freq_band, 'r');

hRev = plot_colored_band(f2(:,1), f2(:,12), freq_band, 'b');
hFwd = plot_colored_band(f2(:,1), f2(:,11), freq_band, 'r');

xlim([0 40]);
ylim([0 max([f2(:,11); f2(:,12)])*1.1]);

xlabel('Frequency (Hz)','FontSize',14);
ylabel('Coherence','FontSize',14);
title('Brain ↔ Spine','FontSize',16);
box off; set(gca,'FontSize',12);
legend([hRev hFwd],{'Reverse','Forward'},'Location','northeast');

% -----------------------
% Spine ↔ EMG
subplot(3,1,3); hold on;

plot_band(f3(:,1), f3(:,12), freq_band, 'b');
plot_band(f3(:,1), f3(:,11), freq_band, 'r');

hRev = plot_colored_band(f3(:,1), f3(:,12), freq_band, 'b');
hFwd = plot_colored_band(f3(:,1), f3(:,11), freq_band, 'r');

xlim([0 40]);
ylim([0 max([f3(:,11); f3(:,12)])*1.1]);

xlabel('Frequency (Hz)','FontSize',14);
ylabel('Coherence','FontSize',14);
title('Spine ↔ EMG','FontSize',16);
box off; set(gca,'FontSize',12);
legend([hRev hFwd],{'Reverse','Forward'},'Location','northeast');



freq_band = [10 35];  % Hz

%% Function to compute forward/reverse area in a given frequency band
compute_area = @(fmat, fwd_col, rev_col, band) ...
    deal( ...
        trapz(fmat(fmat(:,1)>=band(1) & fmat(:,1)<=band(2),1), fmat(fmat(:,1)>=band(1) & fmat(:,1)<=band(2),fwd_col)), ...
        trapz(fmat(fmat(:,1)>=band(1) & fmat(:,1)<=band(2),1), fmat(fmat(:,1)>=band(1) & fmat(:,1)<=band(2),rev_col)) ...
    );

[brainEMG_forward_area, brainEMG_reverse_area] = compute_area(f1, 11, 12, freq_band);
[brainSpine_forward_area, brainSpine_reverse_area] = compute_area(f2, 11, 12, freq_band);
[spineEMG_forward_area, spineEMG_reverse_area] = compute_area(f3, 11, 12, freq_band);


%% Compute directionality ratios
brainEMG_ratio = brainEMG_forward_area / (brainEMG_forward_area + brainEMG_reverse_area);
brainSpine_ratio = brainSpine_forward_area / (brainSpine_forward_area + brainSpine_reverse_area);
spineEMG_ratio = spineEMG_forward_area / (spineEMG_forward_area + spineEMG_reverse_area);


 % Save to structure
    results(ss).sub = sub;

    results(ss).brainEMG.forward_area = brainEMG_forward_area;
    results(ss).brainEMG.reverse_area = brainEMG_reverse_area;
    results(ss).brainEMG.ratio = brainEMG_ratio;

    results(ss).brainSpine.forward_area = brainSpine_forward_area;
    results(ss).brainSpine.reverse_area = brainSpine_reverse_area;
    results(ss).brainSpine.ratio = brainSpine_ratio;

    results(ss).spineEMG.forward_area = spineEMG_forward_area;
    results(ss).spineEMG.reverse_area = spineEMG_reverse_area;
    results(ss).spineEMG.ratio = spineEMG_ratio;


lat_min = -50;
lat_max = 50;

% Brain→EMG
[results(ss).brainEMG.peak_rho, results(ss).brainEMG.peak_latency] = select_peak(t1, lat_min, lat_max);

% Brain→Spine
[results(ss).brainSpine.peak_rho, results(ss).brainSpine.peak_latency] = select_peak(t2, lat_min, lat_max);

% Spine→EMG
[results(ss).spineEMG.peak_rho, results(ss).spineEMG.peak_latency] = select_peak(t3, lat_min, lat_max);

% subplot(3,1,1); hold on
% plot(results(ss).brainEMG.lag, results(ss).brainEMG.rho_time, 'r', 'LineWidth',2);
% xlabel('Lag (ms)'); ylabel('\rho'); title(sprintf('Brain ↔ EMG - %s', sub)); box off; set(gca,'FontSize',12);
end


nSub = length(results);
brainEMG_lat = NaN(nSub,1);
brainSpine_lat = NaN(nSub,1);
spineEMG_lat = NaN(nSub,1);

figure; hold on;

brainEMG_lat   = arrayfun(@(s) s.brainEMG.peak_latency, results);
brainSpine_lat = arrayfun(@(s) s.brainSpine.peak_latency, results);
spineEMG_lat   = arrayfun(@(s) s.spineEMG.peak_latency, results);

% Scatter for individual subjects
scatter(ones(nSub,1), brainEMG_lat, 100, 'filled', 'MarkerFaceAlpha',0.6); % x=1
scatter(2*ones(nSub,1), brainSpine_lat, 100, 'filled', 'MarkerFaceAlpha',0.6); % x=2
scatter(3*ones(nSub,1), spineEMG_lat, 100, 'filled', 'MarkerFaceAlpha',0.6); % x=3

% Compute means and SDs
mean_vals = [mean(brainEMG_lat,'omitnan'), mean(brainSpine_lat,'omitnan'), mean(spineEMG_lat,'omitnan')];
std_vals  = [std(brainEMG_lat,'omitnan'),  std(brainSpine_lat,'omitnan'),  std(spineEMG_lat,'omitnan')];

% Plot mean as large red diamonds
scatter([1 2 3], mean_vals, 200, 'd', 'MarkerFaceColor','k', 'MarkerEdgeColor','k', 'LineWidth',1);

% Add error bars for SD
errorbar([1 2 3], mean_vals, std_vals, 'k','LineWidth',1,'CapSize',10);

% Axis formatting
xlim([0.5 3.5]);
xticks([1 2 3]);
xticklabels({'Brain↔EMG','Brain↔Spine','Spine↔EMG'});
ylabel('Peak latency (ms)');
grid on;
box on;
    ax = gca;              
    ax.FontSize = 14;

%% coherence stuff



nSubs = length(results);

% Preallocate
brainEMG_ratio_all = zeros(nSubs,1);
brainSpine_ratio_all = zeros(nSubs,1);
spineEMG_ratio_all = zeros(nSubs,1);

for ss=1:nSubs
    brainEMG_ratio_all(ss) = results(ss).brainEMG.ratio;
    brainSpine_ratio_all(ss) = results(ss).brainSpine.ratio;
    spineEMG_ratio_all(ss) = results(ss).spineEMG.ratio;
end

% Mean ± SEM
mean_brainEMG = mean(brainEMG_ratio_all);
sem_brainEMG = std(brainEMG_ratio_all)/sqrt(nSubs);

mean_brainSpine = mean(brainSpine_ratio_all);
sem_brainSpine = std(brainSpine_ratio_all)/sqrt(nSubs);

mean_spineEMG = mean(spineEMG_ratio_all);
sem_spineEMG = std(spineEMG_ratio_all)/sqrt(nSubs);

fprintf('Brain→EMG ratio: %.3f ± %.3f\n', mean_brainEMG, sem_brainEMG);
fprintf('Brain→Spine ratio: %.3f ± %.3f\n', mean_brainSpine, sem_brainSpine);
fprintf('Spine→EMG ratio: %.3f ± %.3f\n', mean_spineEMG, sem_spineEMG);

%% Number of subjects
nSubs = numel(results);

%% Initialize matrices
forward_areas = zeros(nSubs, 3); % columns: brain→EMG, brain→spine, spine→EMG
reverse_areas = zeros(nSubs, 3);

for ss = 1:nSubs
    forward_areas(ss,1) = results(ss).brainEMG.forward_area;
    forward_areas(ss,2) = results(ss).brainSpine.forward_area;
    forward_areas(ss,3) = results(ss).spineEMG.forward_area;

    reverse_areas(ss,1) = results(ss).brainEMG.reverse_area;
    reverse_areas(ss,2) = results(ss).brainSpine.reverse_area;
    reverse_areas(ss,3) = results(ss).spineEMG.reverse_area;
end

%% Compute mean and SEM
mean_forward = mean(forward_areas,1);
mean_reverse = mean(reverse_areas,1);
sem_forward = std(forward_areas,0,1)/sqrt(nSubs);
sem_reverse = std(reverse_areas,0,1)/sqrt(nSubs);

%% Prepare data for grouped bar
data = [mean_forward; mean_reverse]'; % rows = connections, columns = Forward/Reverse
sems = [sem_forward; sem_reverse]';

x_labels = {'Brain→EMG','Brain→Spine','Spine→EMG'};

%% Plot grouped bars
figure; hold on;
b = bar(data, 'grouped');
b(1).FaceColor = [0.3 0.6 0.9]; % Forward
b(2).FaceColor = [0.9 0.3 0.3]; % Reverse

% Overlay error bars
[ngroups, nbars] = size(data);
groupwidth = min(0.8, nbars/(nbars+1.5));
for i = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, data(:,i), sems(:,i), 'k.', 'LineWidth',1.5);
end

% Overlay individual subject points
for i = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    for j = 1:ngroups
        if i == 1
            y = forward_areas(:,j);  % Forward
        else
            y = reverse_areas(:,j);  % Reverse
        end
        scatter(repmat(x(j), nSubs,1), y, 50, 'k', 'filled', 'MarkerFaceAlpha',0.6);
    end
end

% Labels and aesthetics
set(gca, 'XTick', 1:ngroups, 'XTickLabel', x_labels, 'FontSize', 14);
ylabel('Coherence area (10-35 Hz)');
ylim([0 max(data(:))*1.3]);
box off;
legend({'Forward','Reverse'}, 'Location', 'northwest');


%% time domain
brainEMG_peaks = zeros(nSubs,1);
brainSpine_peaks = zeros(nSubs,1);
spineEMG_peaks = zeros(nSubs,1);

for ss = 1:nSubs
    brainEMG_peaks(ss) = results(ss).brainEMG.peak_latency;
    brainSpine_peaks(ss) = results(ss).brainSpine.peak_latency;
    spineEMG_peaks(ss) = results(ss).spineEMG.peak_latency;
end

fprintf('Mean ± SEM peak latency (ms):\n');
fprintf('Brain→EMG: %.1f ± %.1f\n', mean(brainEMG_peaks), std(brainEMG_peaks)/sqrt(nSubs));
fprintf('Brain→Spine: %.1f ± %.1f\n', mean(brainSpine_peaks), std(brainSpine_peaks)/sqrt(nSubs));
fprintf('Spine→EMG: %.1f ± %.1f\n', mean(spineEMG_peaks), std(spineEMG_peaks)/sqrt(nSubs));

figure('Color','w'); hold on;

x_labels = {'Brain→EMG','Brain→Spine','Spine→EMG'};
x = 1:3;

% Individual points
scatter(repmat(x(1), nSubs,1), brainEMG_peaks, 50, 'k', 'filled', 'MarkerFaceAlpha',0.6);
scatter(repmat(x(2), nSubs,1), brainSpine_peaks, 50, 'k', 'filled', 'MarkerFaceAlpha',0.6);
scatter(repmat(x(3), nSubs,1), spineEMG_peaks, 50, 'k', 'filled', 'MarkerFaceAlpha',0.6);

% Overlay mean ± SEM
means = [mean(brainEMG_peaks), mean(brainSpine_peaks), mean(spineEMG_peaks)];
sems  = [std(brainEMG_peaks)/sqrt(nSubs), std(brainSpine_peaks)/sqrt(nSubs), std(spineEMG_peaks)/sqrt(nSubs)];
errorbar(x, means, sems, 'r', 'LineWidth',2);

xlim([0.5 3.5]);
set(gca,'XTick',x,'XTickLabel',x_labels,'FontSize',14);
ylabel('Peak latency (ms)');
title('Time-domain directionality peak latencies');
box off;
