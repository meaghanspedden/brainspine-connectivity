%% =========================================================================
% EARLY VS LATE PEAK CMC (first vs last 5 one-second contraction trials)
% =========================================================================

clear all; close all; clc;

%% CONFIG
% Machine-specific paths live in source/brainspine_config.m — edit that
% file to match your local installation.
repo_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(repo_root, 'source'));
paths = brainspine_config();

fieldtrip_path = paths.fieldtrip_path;
spm_path       = paths.spm_path;
neurospec_path = paths.neurospec_path;
data_root      = paths.data_root;
save_dir       = paths.save_dir;

brain_ve_suffix  = '_brain_pct10';

subs_brain = {'OP00212','OP00213','OP00215','OP00219', ...
              'OP00225','OP00221','OP00224'};
subs_spine = {'OP00212','OP00213','OP00215','OP00219', ...
              'OP00220','OP00221','OP00224','OP00225','OP00226'};

subs = intersect(subs_brain,subs_spine,'stable');

fband      = [10 35];
win_sec    = 5;    % number of 1-s contraction trials taken from start/end
seg_pwr_early = 9; % SHORTER than the main analysis (2^11): need enough
                    % complete segments to fit in a 5-s window. 
opt_str    = 'M4';

%% SETUP
addpath(spm_path); spm('defaults','EEG');
addpath(fieldtrip_path); ft_defaults;
addpath(genpath(neurospec_path));

%% STORAGE
peakEarly = nan(numel(subs),1);
peakLate  = nan(numel(subs),1);

%% SUBJECT LOOP
for ss = 1:numel(subs)

    sub = subs{ss};
    fprintf('Processing %s (%d/%d)\n', sub, ss, numel(subs));

    run = '001';
    if strcmp(sub,'OP00224')
        run = '002';
    end

    datafile = fullfile(data_root,['sub-' sub], ...
        'ses-001','meg', ...
        sprintf('pmergedoe1000mspddfflo45hi45hfcstatic_%s_array1.mat',run));

    D = spm_eeg_load(datafile);
    ftdat = spm2fieldtrip(D);
    samp_rate = ftdat.hdr.Fs;

    %% Load VE
    bVE_file = fullfile(save_dir,...
        sprintf('sub%s_VE_brain_M1%s.mat',sub,brain_ve_suffix));
    tmp = load(bVE_file,'VE_brain');
    bVE = tmp.VE_brain;

    %% EMG
    cfg = []; cfg.channel = 'EXG1';
    EMG = ft_selectdata(cfg,ftdat);
    cfg = []; cfg.rectify = 'yes';
    EMG = ft_preprocessing(cfg,EMG);

    %% Contraction indices: first and last win_sec trials
    statidx = find(ftdat.trialinfo == 1);
    early_idx = statidx(1:win_sec);
    late_idx  = statidx(end-win_sec+1:end);

    concat_trials = @(X) cat(2, X.trial{:});

    %% EARLY
    cfg = []; cfg.trials = early_idx;
    B_early   = ft_selectdata(cfg,bVE);
    EMG_early = ft_selectdata(cfg,EMG);
    B_early_cont   = concat_trials(B_early);
    EMG_early_cont = abs(concat_trials(EMG_early));

    %% LATE
    cfg = []; cfg.trials = late_idx;
    B_late   = ft_selectdata(cfg,bVE);
    EMG_late = ft_selectdata(cfg,EMG);
    B_late_cont   = concat_trials(B_late);
    EMG_late_cont = abs(concat_trials(EMG_late));

    %% Check segment count once (first subject only)
    if ss == 1
        n_segs_check = floor(size(B_early_cont,2) / 2^seg_pwr_early);
        fprintf('  seg_pwr_early=%d -> %d samples/segment, %d complete segments in %ds window\n', ...
            seg_pwr_early, 2^seg_pwr_early, n_segs_check, win_sec);
        if n_segs_check < 4
            warning('Fewer than 4 segments — consider lowering seg_pwr_early further.');
        end
    end

    %% Coherence (same estimator as main analysis, shorter segment length)
    [fE, ~, ~] = sp2a2_R2_mt(B_early_cont', EMG_early_cont', samp_rate, seg_pwr_early, opt_str);
    [fL, ~, ~] = sp2a2_R2_mt(B_late_cont',  EMG_late_cont',  samp_rate, seg_pwr_early, opt_str);

    %% Peak coherence (10-35 Hz)
    freq = fE(:,1);
    mask = freq >= fband(1) & freq <= fband(2);

    peakEarly(ss) = max(fE(mask,4));
    peakLate(ss)  = max(fL(mask,4));

end

%% =========================================================================
% STATS
% =========================================================================
[p,~,stats] = signrank(peakEarly, peakLate);

fprintf('\nPeak coherence early vs late:\n');
fprintf('signedrank = %.1f, p = %.4f\n', stats.signedrank, p);

%% =========================================================================
% FIGURE
% =========================================================================
pt_cmap = lines(numel(subs));
p1_idx  = 1;

jit = 0.06;
hfig_early_late = figure('Color','w','Position',[100 100 360 420]); hold on;

for ss = 1:numel(subs)
    xj = [1+jit*randn, 2+jit*randn];
    yj = [peakEarly(ss), peakLate(ss)];

    if ss == p1_idx
        plot(xj, yj, '-o', 'Color','k', 'LineWidth',0.8, 'MarkerSize',10, ...
            'MarkerFaceColor',pt_cmap(ss,:), 'MarkerEdgeColor','k');
    else
        plot(xj, yj, '-', 'Color',[0.75 0.75 0.75], 'LineWidth',0.8);
        scatter(xj(1), yj(1), 55, pt_cmap(ss,:), 'filled', ...
            'MarkerFaceAlpha',0.6, 'MarkerEdgeColor','none');
        scatter(xj(2), yj(2), 55, pt_cmap(ss,:), 'filled', ...
            'MarkerFaceAlpha',0.6, 'MarkerEdgeColor','none');
    end

    text(xj(1)-0.18, yj(1), sprintf('P%d', ss), 'FontSize', 9, ...
        'Color', pt_cmap(ss,:)*0.75, 'HorizontalAlignment','right');
end

med_early = median(peakEarly); med_late = median(peakLate);
mad_early = mad(peakEarly,1);  mad_late  = mad(peakLate,1);
plot(1, med_early, 's','MarkerSize',13,'MarkerFaceColor','w','MarkerEdgeColor','k','LineWidth',1.5);
plot(2, med_late,  's','MarkerSize',13,'MarkerFaceColor','w','MarkerEdgeColor','k','LineWidth',1.5);
errorbar([1 2], [med_early med_late], [mad_early mad_late], 'k','LineWidth',1.5,'CapSize',8,'LineStyle','none');

set(gca,'XTick',[1 2],'XTickLabel',{'Early','Late'},'FontSize',13);
ylabel('Peak Brain-EMG coherence (10-35 Hz)','FontSize',12);
xlim([0.5 2.5]); box off;

text(-0.02, 1.08, 'A', 'Units','normalized', 'FontSize',18, 'FontWeight','bold', ...
    'HorizontalAlignment','left', 'VerticalAlignment','bottom');

fig_dir = fullfile(save_dir, 'figures');
if ~exist(fig_dir,'dir'), mkdir(fig_dir); end
print(hfig_early_late, fullfile(fig_dir,'suppfig_early_vs_late_coherence.png'), '-dpng', '-r300');