%% brainspine_coherence.m
%
%  Full coherence (col 4) for brain-cord, same VEs/settings as
%  RUN_SPECTRA_BS.m. Per-subject spectra are greyed out outside the
%  10-35 Hz band of interest; black stars mark in-band frequencies
%  where full coherence exceeds each subject's own NeuroSpec 95%
%  confidence limit.
%
%  Slot 1 holds a group-level summary panel: each subject's in-band
%  peak coherence (dot) plotted against their own NeuroSpec threshold
%  (black tick), connected by a line and colored by whether the peak
%  clears the bar. This is the per-subject "evidence panel" for the
%  manuscript's Brain-Spine figure; the anatomical ROI recap panel is
%  built separately once a combined brain+spine mesh rendering exists.
%
%  2-row grid layout for paper figure.

clear all; close all; clc;

%% CONFIG (matches RUN_SPECTRA_BS.m)
fieldtrip_path   = 'C:\Users\mspedden\Documents\fieldtrip';
spm_path         = 'C:\Users\mspedden\Documents\spm';
neurospec_path   = 'C:\Users\mspedden\Documents\neurospec211NEW\neurospec211';
data_root        = 'C:\spinecoh_data';
save_dir         = 'C:\Users\mspedden\Documents\brainspine_savetest';

brain_ve_suffix  = '_brain_pct10';
spine_ve_pattern = 'VE_spine_prevalence_sub%s_forspectra_BS.mat';
spine_ve_varname = 'VE';

subs    = {'OP00212','OP00213','OP00215','OP00219', ...
           'OP00225','OP00221','OP00224'};
fband   = [10 35];
seg_pwr = 11;
opt_str = 'M4';

ymax_fixed = 8;                  % Fixed y-axis ceiling (x1e-3 units)
star_y     = 0.92 * ymax_fixed;  % Positioned just below the top border

%% SETUP
addpath(spm_path); spm('defaults','EEG');
addpath(fieldtrip_path); ft_defaults;
addpath(genpath(neurospec_path));

%% SUBJECT LOOP
nSubs = numel(subs);
nrows = 2;
ncols = ceil((nSubs + 1) / nrows); % Total columns required, including the summary panel
hfig  = figure('Color','w','Position',[100 100 ncols*240 nrows*260]);

col_full = [0.8 0.2 0.4];
col_grey = [0.75 0.75 0.75];
col_pass = col_full;   % peak clears subject's own threshold
col_fail = col_grey;   % peak does not clear threshold

peak_val   = nan(1,nSubs);  % in-band peak coherence, per subject
thresh_val = nan(1,nSubs);  % subject's own NeuroSpec 95% CL, per subject

for ss = 1:nSubs
    sub = subs{ss};
    fprintf('Subject %s (%d/%d)\n', sub, ss, nSubs);

    run = '001'; if strcmp(sub,'OP00224'), run = '002'; end
    datafile = fullfile(data_root, ['sub-' sub], 'ses-001', 'meg', ...
        sprintf('pmergedoe1000mspddfflo45hi45hfcstatic_%s_array1.mat', run));
    D     = spm_eeg_load(datafile);
    ftdat = spm2fieldtrip(D);

    bVE_file = fullfile(save_dir, sprintf('sub%s_VE_brain_M1%s.mat', sub, brain_ve_suffix));
    bVE = load(bVE_file, 'VE_brain'); bVE = bVE.VE_brain;

    sVE_file = fullfile(save_dir, sprintf(spine_ve_pattern, sub));
    sVE_data = load(sVE_file); sVE = sVE_data.(spine_ve_varname);

    statidx = find(ftdat.trialinfo==1);
    cfg_ft  = []; cfg_ft.trials = statidx;
    statB   = ft_selectdata(cfg_ft, bVE);
    statS   = ft_selectdata(cfg_ft, sVE);

    statBcont = [statB.trial{:}];
    statScont = [statS.trial{:}];
    samp_rate = ftdat.hdr.Fs;

    [f_mt, ~, cl_mt] = sp2a2_R2_mt(statBcont', statScont', samp_rate, seg_pwr, opt_str);

    freq_axis = f_mt(:,1)';
    coh_full  = f_mt(:,4)' * 1e3;
    sig_line  = cl_mt.ch_c95 * 1e3;

    mask_in_band  = freq_axis >= fband(1) & freq_axis <= fband(2);
    mask_out_band = ~mask_in_band;

    coh_full_grey  = coh_full; coh_full_grey(~mask_out_band) = NaN;
    coh_full_color = coh_full; coh_full_color(~mask_in_band)  = NaN;

    % In-band peak coherence and this subject's own significance bar
    peak_val(ss)   = max(coh_full(mask_in_band));
    thresh_val(ss) = sig_line;

    % Significance restricted to the 10-35 Hz band of interest
    sig_any = (coh_full > sig_line) & mask_in_band;

    % Plotting index shifted by +1 to leave subplot(nrows, ncols, 1) for the summary panel
    plot_idx = ss + 1;
    subplot(nrows, ncols, plot_idx); hold on;

    % ROI Shading
    fill([fband(1) fband(2) fband(2) fband(1)], [0 0 ymax_fixed ymax_fixed], ...
        [0.95 0.95 0.95], 'FaceAlpha',0.5, 'EdgeColor','none', 'HandleVisibility','off');

    % Data Plots
    plot(freq_axis, coh_full_grey,  '-', 'Color',col_grey, 'LineWidth',1,    'HandleVisibility','off');
    plot(freq_axis, coh_full_color, '-', 'Color',col_full, 'LineWidth',1.75, 'HandleVisibility','off');

    % Significance Stars
    if any(sig_any)
        plot(freq_axis(sig_any), star_y*ones(1,sum(sig_any)), '*', ...
            'Color','k', 'MarkerSize',4, 'HandleVisibility','off');
    end

    % Fine-Tuning Axes Aesthetics
    xlim([2 45]); ylim([0 ymax_fixed]);
    set(gca, 'XTick', 10:10:40, 'YTick', 0:2:ymax_fixed);
    set(gca, 'FontSize', 9, 'FontName', 'Helvetica', 'TickDir', 'out', 'LineWidth', 0.8);
    box off;

    if mod(plot_idx-1, ncols) == 0
        ylabel('Coherence (\times10^{-3})', 'FontSize', 10, 'FontWeight', 'bold');
    else
        set(gca, 'YTickLabel', []);
    end

    if plot_idx > (nrows-1)*ncols
        xlabel('Frequency (Hz)', 'FontSize', 10, 'FontWeight', 'bold');
    else
        set(gca, 'XTickLabel', []);
    end

    title(sprintf('Participant %d', ss), 'FontSize', 10, 'FontWeight', 'bold', 'Color', [0.2 0.2 0.2]);
end

%% GROUP SUMMARY PANEL (slot 1): peak coherence vs each subject's own threshold
pass_flag = peak_val > thresh_val;
tick_hw   = 0.18; % half-width of the threshold tick mark, in x-axis units

subplot(nrows, ncols, 1); hold on;

for ss = 1:nSubs
    if pass_flag(ss), c = col_pass; else, c = col_fail; end
    plot([ss ss], [thresh_val(ss) peak_val(ss)], '-', 'Color', c, 'LineWidth', 1.5, 'HandleVisibility','off');
    plot([ss-tick_hw ss+tick_hw], [thresh_val(ss) thresh_val(ss)], '-', 'Color','k', 'LineWidth',2, 'HandleVisibility','off');
    plot(ss, peak_val(ss), 'o', 'MarkerFaceColor', c, 'MarkerEdgeColor','k', 'MarkerSize',7, 'HandleVisibility','off');
end

% Dummy handles, for a compact legend
h_pass = plot(nan, nan, 'o', 'MarkerFaceColor',col_pass, 'MarkerEdgeColor','k', 'MarkerSize',7);
h_fail = plot(nan, nan, 'o', 'MarkerFaceColor',col_fail, 'MarkerEdgeColor','k', 'MarkerSize',7);
h_thr  = plot([NaN NaN], [NaN NaN], '-', 'Color','k', 'LineWidth',2);
legend([h_pass h_fail h_thr], {'Clears threshold','Below threshold','Subject threshold'}, ...
    'Location','best', 'FontSize',7, 'Box','off');

xlim([0.5 nSubs+0.5]); ylim([0 ymax_fixed]);
set(gca, 'XTick', 1:nSubs, 'XTickLabel', arrayfun(@(x) sprintf('P%d',x), 1:nSubs, 'UniformOutput', false));
set(gca, 'YTick', 0:2:ymax_fixed);
set(gca, 'FontSize', 9, 'FontName', 'Helvetica', 'TickDir', 'out', 'LineWidth', 0.8);
box off;

ylabel('Coherence (\times10^{-3})', 'FontSize', 10, 'FontWeight', 'bold');
xlabel('Participant', 'FontSize', 10, 'FontWeight', 'bold');
title('Peak Coherence vs. Threshold', 'FontSize', 10, 'FontWeight', 'bold', 'Color', [0.2 0.2 0.2]);

sgtitle('Brain-Cord Coherence Profiles (* = p < 0.05 vs NeuroSpec, 10-35 Hz)', ...
        'FontSize', 12, 'FontWeight', 'bold', 'FontName', 'Helvetica');

fig_dir = fullfile(save_dir, 'figures');
if ~exist(fig_dir,'dir'), mkdir(fig_dir); end

set(hfig, 'Units', 'Inches');
pos = get(hfig, 'Position');
set(hfig, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Inches', 'PaperSize', [pos(3), pos(4)]);

print(hfig, fullfile(fig_dir, 'brainspine_coherence_spectra_and_summary.png'), '-dpng', '-r300');
fprintf('Figure saved successfully.\n');