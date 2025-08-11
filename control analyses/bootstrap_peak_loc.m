% Bootstrap analysis of spinal source peak power 

close all
clear all
clc

addpath('D:\spm')
spm('defaults','EEG')

subs = {'OP00212', 'OP00213',  'OP00215', 'OP00219', ...
   'OP00220', 'OP00221',  'OP00225', 'OP00226'};

%subs={'OP00214'}; %d
%subs={'OP00224'};%002
freqband=[10 35];
which_ori='all'; % 'all', '1', '2', '3'
n_bootstrap = 500; % Number of bootstrap iterations

addpath('D:\brainspineconnectivity\stats')

%% filenames for brain spinal cord and emg------------------------------------------
for s=1:length(subs)
    sub=subs{s};
save_dir = fullfile('D:\MSST001', [sub '_contrast']);
if ~exist(save_dir, 'dir')
    mkdir(save_dir)
end
bffilename = fullfile(save_dir, ['bfdata_', sub]);

datwithEMGmerged = fullfile('D:\MSST001', ...
    ['sub-' sub], ...
    'ses-001', ...
    'meg', ...
    ['pmergedoe1000mspddfflo45hi45hfcstatic_001_array1.mat']);

gemfile = 'D:\MSST001\generic_merged\geoms.mat';

%% load data
D=spm_eeg_load(datwithEMGmerged);
grad=D.sensors('MEG');

castforward = load(gemfile);
mtorso = castforward.mesh_torso;
src=castforward.sources_center_line;

st=load(bffilename);

dat=spm2fieldtrip(D);
num_sources=size(st.wdatastatic,1);

cfg=[];
cfg.channel = dat.label(1:num_sources);
alldat=ft_selectdata(cfg,dat);

cfg=[];
cfg.channel='EXG1';
cfg.rectify ='yes';
cfg.detrend='yes';
EMG=ft_preprocessing(cfg, dat);

alldat_EMG=ft_appenddata([], alldat, EMG);

if strcmp(which_ori, 'all')
    ori_list = 1:3;
else
    ori_list = str2double(which_ori);
end

stattrials=find(strcmp(D.conditions,'static'));
resttrials=find(strcmp(D.conditions,'rest'));
    
ntrials = min([numel(stattrials) numel(resttrials)]);

all_powstat = zeros(ntrials, num_sources);
all_powrest = zeros(ntrials, num_sources);

for o = ori_list
    for k = 1:numel(alldat.trial)
        alldat.trial{k}(:,:) = squeeze(st.wdata(:,k,:,o));
    end

    if o == ori_list(1)
        for j = 1:numel(alldat.label)-1
            alldat.label{j} = sprintf('source%g', j);
        end
    end

    cfg = [];
    cfg.output     = 'pow';
    cfg.method     = 'mtmfft';
    cfg.foilim     = freqband;
    cfg.tapsmofrq  = 2;
    cfg.keeptrials = 'yes';
    cfg.pad        = 'nextpow2';

    fdat = ft_freqanalysis(cfg, alldat);

    meanpowstat = mean(fdat.powspctrm(stattrials,:,:), 3);
    meanpowrest = mean(fdat.powspctrm(resttrials,:,:), 3);

    if strcmp(which_ori, 'all')
        all_powstat = all_powstat + meanpowstat(1:ntrials,:);
        all_powrest = all_powrest + meanpowrest(1:ntrials,:);
    else
        all_powstat = meanpowstat;
        all_powrest = meanpowrest;
    end
end

%% BOOTSTRAP PEAK LOCATION 
peak_indices = zeros(1, n_bootstrap);
bootstrap_peak_tvals = nan(1, n_bootstrap);

diffmat = log(all_powstat(1:ntrials,:)) - log(all_powrest(1:ntrials,:));

for i = 1:n_bootstrap
    % Resample trials with replacement
    resample_idx = randsample(ntrials, ntrials, true);
    
    % Compute average difference across resampled trials
    resampled_diff = diffmat(resample_idx,:);
    avg_diff = mean(resampled_diff, 1);
    
    % Find peak index (min avg_diff)
    [~, peak_idx] = min(avg_diff);
    peak_indices(i) = peak_idx;
    
    % Get data at bootstrap peak
    stat_vals = log(all_powstat(resample_idx, peak_idx));
    rest_vals = log(all_powrest(resample_idx, peak_idx));
    
    % Perform one-tailed t-test (left tail)
    [~, ~, ~, stats] = ttest2(stat_vals, rest_vals, 'Tail', 'left');
    bootstrap_peak_tvals(i) = stats.tstat;
end

% Original peak
orig_diff = mean(log(all_powstat)) - mean(log(all_powrest));
[~, orig_peak_idx] = min(orig_diff);
orig_peak_y = src.pos(orig_peak_idx, 2);

% Compute original peak t-stat
[~, ~, ~, stats] = ttest2(log(all_powstat(:, orig_peak_idx)), log(all_powrest(:, orig_peak_idx)), 'Tail', 'left');
orig_t = stats.tstat;


%% PLOT PEAK HISTOGRAM

% Get Y positions of all bootstrap peaks
% Bootstrap Y positions (example data)
bootstrap_y_positions = src.pos(peak_indices, 2);

% Median and MAD
med_y = median(bootstrap_y_positions);
mad_y = mad(bootstrap_y_positions, 1);
mad_std = max(mad_y * 1.4826, 1e-2);  % avoid too narrow Gaussian

% Histogram parameters
nBins = 30;
[counts, edges] = histcounts(bootstrap_y_positions, nBins);
bin_centers = (edges(1:end-1) + edges(2:end)) / 2;
bin_width = edges(2) - edges(1);

% Plot histogram
figure;
yyaxis left
bar(bin_centers, counts, 'FaceAlpha', 0.5);
ylabel('Count');
xlabel('Y Position');

% Gaussian x range focused around median ±4 std
% x_gauss = linspace(med_y - 4*mad_std, med_y + 4*mad_std, 200);

% % Gaussian PDF
% gauss_pdf = (1 / (mad_std * sqrt(2*pi))) * exp(-0.5 * ((x_gauss - med_y)/mad_std).^2);
% 
% % Scale Gaussian to histogram counts
% gauss_pdf_scaled = gauss_pdf * length(bootstrap_y_positions) * bin_width;
% 
% % Plot Gaussian on right axis for better visualization
% yyaxis right
% plot(x_gauss, gauss_pdf_scaled, 'r-', 'LineWidth', 2);
% ylabel('Scaled Gaussian');

% Vertical median line
hold on;
xline(med_y, 'k--', 'LineWidth', 2, 'Label', 'Median', 'LabelHorizontalAlignment', 'left');

title(sprintf('Histogram sub %s',sub));
legend({'Histogram', 'Gaussian fit', 'Median'});



% Save if desired
savename = sprintf('spine_bootstrap_%s_%s', sub, which_ori);
save(fullfile(save_dir,savename),  'bootstrap_y_positions', 'med_y','mad_std' )
end