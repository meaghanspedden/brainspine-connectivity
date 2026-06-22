%% brainspine_peak_vs_threshold_boxplot.m
%
%  Standalone figure: a single boxplot (median, IQR, whiskers to
%  min/max) summarizing the distribution of each subject's peak
%  full coherence normalized to their own NeuroSpec 95% confidence
%  limit (peak / threshold). Individual subjects are overlaid as
%  jittered dots, colored by whether they clear the bar (ratio > 1).
%  A faint dashed line at ratio = 1 marks the shared threshold.
%
%  Self-contained: recomputes brain-cord coherence per subject rather
%  than depending on output from another script.

clear all; close all; clc;

%% CONFIG (matches RUN_SPECTRA_BS.m / brainspine_coherence_lag_compare.m)
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

%% SETUP
addpath(spm_path); spm('defaults','EEG');
addpath(fieldtrip_path); ft_defaults;
addpath(genpath(neurospec_path));

%% SUBJECT LOOP: extract peak coherence and threshold, no spectra plotting here
nSubs = numel(subs);
peak_val   = nan(1,nSubs);
thresh_val = nan(1,nSubs);

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

    mask_in_band = freq_axis >= fband(1) & freq_axis <= fband(2);

    peak_val(ss)   = max(coh_full(mask_in_band));
    thresh_val(ss) = sig_line;
end

ratio     = peak_val ./ thresh_val;   % peak coherence normalized to each subject's own threshold
pass_flag = ratio > 1;

%% FIGURE: boxplot of the ratio, individual subjects overlaid and colored by significance
col_pass = [0.8 0.2 0.4];     % matches col_full in the spectra figure
col_fail = [0.65 0.65 0.65];
col_edge = [0.15 0.15 0.15];
col_box  = [0.93 0.93 0.93];
col_line = [0.20 0.20 0.20];
col_thr  = [0.10 0.10 0.10];
col_grid = [0.88 0.88 0.88];

q1  = quantile(ratio, 0.25);
q3  = quantile(ratio, 0.75);
med = median(ratio);
wlo = min(ratio);
whi = max(ratio);

box_hw = 0.22;   % half-width of the box
cap_hw = 0.10;   % half-width of the whisker caps
jitter = linspace(-0.16, 0.16, nSubs);   % deterministic spread, so subjects don't sit in one stack

y_top = max([ratio(:); 1]) * 1.15;
y_top = ceil(y_top*2)/2;   % round up to the nearest 0.5 for clean ticks

hfig = figure('Color','w','Position',[100 100 600 480]);
ax   = axes('Parent', hfig); hold(ax,'on');

% for gy = 0:0.5:y_top
%     plot(ax, [0.3 1.7], [gy gy], '-', 'Color', col_grid, 'LineWidth', 1, 'HandleVisibility','off');
% end

% Shared threshold reference line
plot(ax, [0.3 1.7], [1 1], '--', 'Color', col_thr, 'LineWidth', 1.5, 'HandleVisibility','off');

% Whiskers (min/max, since n is small the 1.5xIQR convention isn't very meaningful here)
plot(ax, [1 1], [q3 whi], '-', 'Color', col_line, 'LineWidth', 1.2, 'HandleVisibility','off');
plot(ax, [1 1], [wlo q1], '-', 'Color', col_line, 'LineWidth', 1.2, 'HandleVisibility','off');
plot(ax, [1-cap_hw 1+cap_hw], [whi whi], '-', 'Color', col_line, 'LineWidth', 1.2, 'HandleVisibility','off');
plot(ax, [1-cap_hw 1+cap_hw], [wlo wlo], '-', 'Color', col_line, 'LineWidth', 1.2, 'HandleVisibility','off');

% Box (Q1-Q3)
fill(ax, [1-box_hw 1+box_hw 1+box_hw 1-box_hw], [q1 q1 q3 q3], col_box, ...
    'EdgeColor', col_line, 'LineWidth', 1.2, 'HandleVisibility','off');

% Median line
plot(ax, [1-box_hw 1+box_hw], [med med], '-', 'Color', col_line, 'LineWidth', 2.2, 'HandleVisibility','off');

% Individual subjects, jittered, colored by significance
for ss = 1:nSubs
    if pass_flag(ss), c = col_pass; else, c = col_fail; end
    plot(ax, 1+jitter(ss), ratio(ss), 'o', 'MarkerFaceColor', c, 'MarkerEdgeColor', col_edge, ...
        'MarkerSize', 9, 'LineWidth', 1, 'HandleVisibility','off');
end

% Highlight a specific participant with an open black ring (e.g. for talk slides)
highlight_idx = 1;   % index into `subs`; 1 = P1
plot(ax, 1+jitter(highlight_idx), ratio(highlight_idx), 'o', 'MarkerFaceColor', 'none', ...
    'MarkerEdgeColor', 'k', 'MarkerSize', 18, 'LineWidth', 2, 'HandleVisibility','off');

% Dummy handles for the legend
h_pass = plot(ax, nan, nan, 'o', 'MarkerFaceColor', col_pass, 'MarkerEdgeColor', col_edge, 'MarkerSize', 9, 'LineWidth', 1);
h_fail = plot(ax, nan, nan, 'o', 'MarkerFaceColor', col_fail, 'MarkerEdgeColor', col_edge, 'MarkerSize', 9, 'LineWidth', 1);
h_thr  = plot(ax, [NaN NaN], [NaN NaN], '--', 'Color', col_thr, 'LineWidth', 1.5);
legend(ax, [h_pass h_fail h_thr], {'Significant ','Nonsignificant','Threshold'}, ...
    'Location','best', 'FontSize', 10, 'Box','off');

xlim(ax, [0.3 1.7]); ylim(ax, [0 y_top]);
set(ax, 'XTick', 1, 'XTickLabel', sprintf('Participants'));
set(ax, 'YTick', 0:0.5:y_top);
set(ax, 'FontSize', 11, 'FontName', 'Helvetica', 'TickDir', 'out', 'LineWidth', 1);
box(ax, 'off');

ylabel(ax, 'Peak coherence / threshold', 'FontSize', 13, 'FontWeight', 'bold');
title(ax, 'Brain-spinal cord coherence', 'FontSize', 15, 'FontWeight', 'bold', 'Color', [0.15 0.15 0.15]);

fig_dir = fullfile(save_dir, 'figures');
if ~exist(fig_dir,'dir'), mkdir(fig_dir); end

set(hfig, 'Units', 'Inches');
pos = get(hfig, 'Position');
set(hfig, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Inches', 'PaperSize', [pos(3), pos(4)]);

print(hfig, fullfile(fig_dir, 'brainspine_peak_vs_threshold_boxplot.png'), '-dpng', '-r300');
fprintf('Figure saved successfully.\n');