%% peak freq correlation plots
%
%  Peak frequency, line of equality, for the three coherence pairs
%  (Brain-EMG, Brain-Cord, Cord-EMG), using full coherence (sp2a2_R2_mt
%  column 4, all lags including zero) within the 10-35 Hz analysis band.
%  Replaces the phase-clustering / per-trial Rayleigh-test machinery from
%  phase_analysis_v2.m, which is dropped entirely.
%
%  For each participant, find each pair's own peak coherence frequency
%  in-band, then test pairwise (BE vs BS, BE vs SE, BS vs SE) whether the
%  peaks coincide across the group via a paired t-test on the differences
%  (H0: mean diff = 0, i.e. points lie on the y=x line).

% =========================================================================

clear all; close all; clc;

%% =========================================================================
%  CONFIG
%% =========================================================================
fieldtrip_path = 'C:\Users\mspedden\Documents\fieldtrip';
spm_path       = 'C:\Users\mspedden\Documents\spm';
neurospec_path = 'C:\Users\mspedden\Documents\neurospec211NEW\neurospec211';
data_root      = 'C:\spinecoh_data';
save_dir       = 'C:\Users\mspedden\Documents\brainspine_savetest';

brain_ve_suffix  = '_brain_pct10';
spine_ve_pattern = 'VE_spine_prevalence_sub%s_forspectra_BS.mat';
spine_ve_varname = 'VE';

subs_brain = {'OP00212','OP00213','OP00215','OP00219', ...
              'OP00225','OP00221','OP00224'};
subs_spine = {'OP00212','OP00213','OP00215','OP00219', ...
              'OP00220','OP00221','OP00224','OP00225','OP00226'};
subs       = intersect(subs_brain, subs_spine, 'stable');

fband   = [10 35];
seg_pwr = 11;
opt_str = 'M4';

col_BE = [0.2 0.4 0.8];
col_BS = [0.8 0.2 0.4];
col_SE = [0.2 0.7 0.3];

fig_dir = fullfile(save_dir, 'figures', 'peak_freq_and_shape_BS');
if ~exist(fig_dir,'dir'), mkdir(fig_dir); end

%% =========================================================================
%  SETUP
%% =========================================================================
addpath(spm_path); spm('defaults','EEG');
addpath(fieldtrip_path); ft_defaults;
addpath(genpath(neurospec_path));

%% =========================================================================
%  SUBJECT LOOP
%% =========================================================================
nSubs   = numel(subs);
results = struct();

for ss = 1:nSubs
    sub = subs{ss};
    fprintf('\n  Subject %s (%d/%d)\n', sub, ss, nSubs);

    %% Load data
    run = '001'; if strcmp(sub,'OP00224'), run = '002'; end
    datafile = fullfile(data_root, ['sub-' sub], 'ses-001', 'meg', ...
        sprintf('pmergedoe1000mspddfflo45hi45hfcstatic_%s_array1.mat', run));
    D         = spm_eeg_load(datafile);
    ftdat     = spm2fieldtrip(D);
    samp_rate = ftdat.hdr.Fs;

    %% Load VEs
    bVE_file = fullfile(save_dir, sprintf('sub%s_VE_brain_M1%s.mat', sub, brain_ve_suffix));
    sVE_file = fullfile(save_dir, sprintf(spine_ve_pattern, sub));
    bVE      = load(bVE_file, 'VE_brain'); bVE = bVE.VE_brain;
    sVE_data = load(sVE_file); sVE = sVE_data.(spine_ve_varname);

    %% EMG
    cfg_ft = []; cfg_ft.channel = 'EXG1';
    EMG = ft_selectdata(cfg_ft, ftdat);
    cfg_ft = []; cfg_ft.rectify = 'yes';
    EMG = ft_preprocessing(cfg_ft, EMG);

    %% Contraction trials
    statidx = find(ftdat.trialinfo == 1);
    cfg_ft  = []; cfg_ft.trials = statidx;
    statB   = ft_selectdata(cfg_ft, bVE);
    statS   = ft_selectdata(cfg_ft, sVE);
    statEMG = ft_selectdata(cfg_ft, EMG);

    statBcont   = [statB.trial{:}];
    statScont   = [statS.trial{:}];
    statEMGcont = abs([statEMG.trial{:}]);

    %% Coherence spectra (multitaper, full coherence = column 4)
    [f_BE, ~, ~] = sp2a2_R2_mt(statBcont', statEMGcont', samp_rate, seg_pwr, opt_str);
    [f_BS, ~, ~] = sp2a2_R2_mt(statBcont', statScont',   samp_rate, seg_pwr, opt_str);
    [f_SE, ~, ~] = sp2a2_R2_mt(statScont', statEMGcont', samp_rate, seg_pwr, opt_str);

    freq_axis  = f_BE(:,1)';
    band_mask  = freq_axis >= fband(1) & freq_axis <= fband(2);
    freqs_band = freq_axis(band_mask);

    coh_BE = f_BE(band_mask, 4);
    coh_BS = f_BS(band_mask, 4);
    coh_SE = f_SE(band_mask, 4);

    %% Peak frequencies (each pair independently, full coherence)
    [~, ipk_BE] = max(coh_BE);
    [~, ipk_BS] = max(coh_BS);
    [~, ipk_SE] = max(coh_SE);

    peak_BE = freqs_band(ipk_BE);
    peak_BS = freqs_band(ipk_BS);
    peak_SE = freqs_band(ipk_SE);

    fprintf('  Peak freqs (full coherence) - BE: %.1f Hz  BS: %.1f Hz  SE: %.1f Hz\n', ...
        peak_BE, peak_BS, peak_SE);

    results(ss).sub        = sub;
    results(ss).peak_BE    = peak_BE;
    results(ss).peak_BS    = peak_BS;
    results(ss).peak_SE    = peak_SE;
    results(ss).coh_BE     = coh_BE;
    results(ss).coh_BS     = coh_BS;
    results(ss).coh_SE     = coh_SE;
    results(ss).freqs_band = freqs_band;
end

%% =========================================================================
%  SAVE
%% =========================================================================
save(fullfile(save_dir, 'peak_freq_and_shape_BS_results.mat'), 'results');

%% =========================================================================
%  SUMMARY + FIGURE - peak frequency, line of equality
%% =========================================================================
peak_BE_all = [results.peak_BE]';
peak_BS_all = [results.peak_BS]';
peak_SE_all = [results.peak_SE]';

comparisons = {
    peak_BE_all, peak_BS_all, 'Brain–EMG peak (Hz)', 'Brain–Cord peak (Hz)', col_BE, col_BS;
    peak_BE_all, peak_SE_all, 'Brain–EMG peak (Hz)', 'Cord–EMG peak (Hz)',   col_BE, col_SE;
    peak_BS_all, peak_SE_all, 'Brain–Cord peak (Hz)', 'Cord–EMG peak (Hz)', col_BS, col_SE};

pt_cmap = lines(nSubs);

hfig1 = figure('Color','w','Position',[100 100 1050 370]);

for pp = 1:3
    xvals   = comparisons{pp,1};
    yvals   = comparisons{pp,2};
    xlbl    = comparisons{pp,3};
    ylbl    = comparisons{pp,4};
    col_fit = (comparisons{pp,5} + comparisons{pp,6}) / 2;

    subplot(1,3,pp); hold on;

    ax_min = min([fband(1); xvals; yvals]) - 1;
    ax_max = max([fband(2); xvals; yvals]) + 1;

    for ss = 1:nSubs
        scatter(xvals(ss), yvals(ss), 90, 'o', ...
            'MarkerFaceColor', pt_cmap(ss,:), ...
            'MarkerEdgeColor', 'w', 'LineWidth', 1.0, ...
            'HandleVisibility', 'off');
        text(xvals(ss) + 0.25, yvals(ss), sprintf('P%d', ss), ...
            'FontSize', 9, 'Color', pt_cmap(ss,:) * 0.75);
    end

    plot([ax_min ax_max], [ax_min ax_max], '--', ...
        'Color', [0.5 0.5 0.5], 'LineWidth', 1.5, 'HandleVisibility', 'off');

    diffs = xvals - yvals;
    [~, t_pval, ~, t_stats] = ttest(diffs);

    stat_str = sprintf('t(%d) = %.2f, p %s', ...
        t_stats.df, t_stats.tstat, pfmt(t_pval));

    text(ax_min + 0.05 * range([ax_min ax_max]), ...
         ax_max - 0.08 * range([ax_min ax_max]), ...
         stat_str, ...
         'FontSize', 10, 'Color', col_fit * 0.8, 'FontWeight', 'bold', ...
         'VerticalAlignment', 'top');

    xlim([ax_min ax_max]); ylim([ax_min ax_max]);
    axis square;
    xlabel(xlbl, 'FontSize', 12);
    ylabel(ylbl, 'FontSize', 12);
    set(gca, 'FontSize', 12, 'GridLineStyle', ':', 'GridAlpha', 0.2);
    box off; grid on;
end

sgtitle('Peak coherence frequency per signal pair (full coherence), all participants', ...
    'FontSize', 13, 'FontWeight', 'normal');
print(hfig1, fullfile(fig_dir, 'group_peak_freq_scatter.png'), '-dpng', '-r300');
fprintf('Peak frequency figure saved.\n');

%% =========================================================================
%  HELPERS
%% =========================================================================
function s = pfmt(p)
    if p < 0.001
        s = '< .001';
    else
        s = sprintf('= %.3f', p);
    end
end