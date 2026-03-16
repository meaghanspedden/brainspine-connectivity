%% RUN_SPECTRA.m
% Downstream spectra and coherence analysis for brain-spine pipeline.
%
% REQUIRES outputs from RUN_PIPELINE.m:
%   - M1_ROI_bemv2_*.mat          (Step 4) — brain VE ROI positions
%   - VE_spine_sub*_bemv2*.mat    (Step 3) — spinal virtual electrodes
%
% STEPS:
%   Step A — Build brain virtual electrode (LCMV) from M1 ROI
%   Step B — Pairwise coherence spectra (NeuroSpec): Brain-EMG, Brain-Spine, Spine-EMG
%             + directed coherence (forward/reverse areas)
%             + FieldTrip coherence spectra
%             + spectral power (10-35 Hz) static vs rest
%
% EDIT USER CONFIG BELOW THEN RUN.
% =========================================================================

clear all
close all
clc

%% =========================================================================
%  USER CONFIG
%  =========================================================================

% --- Toolbox paths --------------------------------------------------------
cfg.fieldtrip_path  = 'C:\Users\mspedden\Documents\fieldtrip';
cfg.spm_path        = 'C:\Users\mspedden\Documents\spm';
cfg.bsc_source_path = 'C:\Users\mspedden\Documents\brainspineconnectivity\source';
cfg.neurospec_path  = 'C:\Users\mspedden\Documents\neurospec211NEW\neurospec211';

% --- Data & geometry paths ------------------------------------------------
cfg.data_root = 'C:\Users\mspedden\Documents';
cfg.geomfile  = 'C:\Users\mspedden\Documents\new_leadfields_and_geom\geometries_cervical_realistic.mat';
cfg.save_dir  = 'C:\Users\mspedden\Documents\brainspine_save_bemv2';

% --- Subject lists --------------------------------------------------------
% Brain VE: 7 subjects (exclude OP00220 small head, OP00226 no headcast)
cfg.subs_brain = {'OP00212','OP00213','OP00215','OP00219', ...
                  'OP00225','OP00221','OP00224'};

% Spine VE: all 9 subjects
cfg.subs_spine = {'OP00212','OP00213','OP00215','OP00219', ...
                  'OP00220','OP00221','OP00224','OP00225','OP00226'};

% Coherence analysis: only subjects with BOTH VEs (intersection)
cfg.subs_coh = intersect(cfg.subs_brain, cfg.subs_spine);

% --- Smoothing option -----------------------------------------------------
% Must match what was used in RUN_PIPELINE to select correct VE files.
% doSmooth = 1 loads VE_spine_sub*_forspectra_bemv2_permSmooth_40mm.mat
% doSmooth = 0 loads VE_spine_sub*_forspectra_bemv2.mat

cfg.doSmooth             = 1;
cfg.spine_smooth_fwhm_mm = 20;   % must match RUN_PIPELINE

% M1 ROI file suffix — must match result_suffix from RUN_PIPELINE Step 4
% e.g. 'functionalVE_spineSmooth_40mm_brainSmooth_8mm'
cfg.M1_roi_suffix = 'functionalVE_spineSmooth_20mm_brainSmooth_8mm';

% --- Analysis parameters --------------------------------------------------
cfg.fband    = [10 35];   % beta band for power/coherence summary
cfg.seg_pwr  = 11;        % NeuroSpec segment power (2^11 = 2048 samples)
cfg.lat_min  = -50;       % cross-correlation latency window (ms)
cfg.lat_max  =  50;

% --- Figure saving -------------------------------------------------------
saveFigs = 1;   % 1 = save all figures as .png to save_dir/figures/

% --- Steps to run ---------------------------------------------------------
run_stepA = 0;   % Brain VE
run_stepB = 1;   % Coherence spectra

%% =========================================================================
%  SETUP
%% =========================================================================

addpath(cfg.bsc_source_path)
addpath(cfg.spm_path)
spm('defaults','EEG')
addpath(cfg.fieldtrip_path)
ft_defaults
addpath(genpath(cfg.neurospec_path))

cfg.saveFigs = saveFigs;

% Build spine VE filename suffix from config
if cfg.doSmooth
    cfg.spine_ve_suffix = sprintf('_permSmooth_%dmm', cfg.spine_smooth_fwhm_mm);
else
    cfg.spine_ve_suffix = '';
end

% Load M1 ROI file
roi_file = fullfile(cfg.save_dir, sprintf('M1_ROI_bemv2_%s.mat', cfg.M1_roi_suffix));
assert(exist(roi_file,'file')==2, ...
    'M1 ROI file not found:\n  %s\nRun RUN_PIPELINE Step 4 first.', roi_file)
load(roi_file)  % loads intersection_pos, roi_native, M1_mni, radius

fprintf('\n=== SPECTRA CONFIG ===\n')
if cfg.doSmooth
    fprintf('  Spine VE:   smoothed (%d mm)\n', cfg.spine_smooth_fwhm_mm)
else
    fprintf('  Spine VE:   unsmoothed\n')
end
fprintf('  Brain subs: %d  |  Spine subs: %d  |  Coh subs: %d\n', ...
    numel(cfg.subs_brain), numel(cfg.subs_spine), numel(cfg.subs_coh))
fprintf('  M1 ROI:     %d positions\n', size(intersection_pos,1))
fprintf('======================\n\n')

%% =========================================================================
%  STEP A — Brain virtual electrode (LCMV)
%% =========================================================================

if run_stepA
    fprintf('\n>>> STEP A: Brain virtual electrode\n')
    run_brain_VE(cfg, intersection_pos)
    fprintf('>>> STEP A complete.\n\n')
end

%% =========================================================================
%  STEP B — Coherence spectra (NeuroSpec + FieldTrip)
%% =========================================================================

if run_stepB
    % check at least first subject has brain VE
    sub1    = cfg.subs_coh{1};
    ve_file = fullfile(cfg.save_dir, sprintf('sub%s_VE_brain_forspectra.mat', sub1));
    assert(exist(ve_file,'file')==2, ...
        'Brain VE not found for %s:\n  %s\nRun Step A first.', sub1, ve_file)
    fprintf('\n>>> STEP B: Coherence spectra\n')
    run_coherence_spectra(cfg)
    fprintf('>>> STEP B complete.\n\n')
end

fprintf('\n=== SPECTRA FINISHED ===\n\n')


%% =========================================================================
%  STEP FUNCTIONS
%% =========================================================================

function run_brain_VE(cfg, intersection_pos)
% Build brain virtual electrode from M1 ROI using LCMV beamformer.

subs     = cfg.subs_brain;
save_dir = cfg.save_dir;
geomfile = cfg.geomfile;
data_root = cfg.data_root;

for ss = 1:length(subs)
    sub = subs{ss};
    fprintf('  [Step A] Subject %s (%d/%d)\n', sub, ss, length(subs));

    datwithEMGmerged = get_datafile(data_root, sub);
    load(geomfile)   % sources_brain, mesh_brain, etc.
    D       = spm_eeg_load(datwithEMGmerged);
    grad_mm = D.sensors('MEG');
    ftdat   = spm2fieldtrip(D);

    badchans = D.chanlabels(D.badchannels);
    cfg_ft = []; cfg_ft.channel = setdiff(ftdat.label, badchans);
    ftdat = ft_selectdata(cfg_ft, ftdat);

    % brain channel subset
    brainidx           = find(grad_mm.chanpos(:,2) > 200);
    braingrad          = subset_grad(grad_mm, brainidx);

    % braindat = full data (matches original VE_brain_spectra.m)
    braindat        = ftdat;
    mesh_brain.unit = 'mm';

    % volume conductor + leadfield
    cfg_ft = []; cfg_ft.method = 'singleshell';
    vol = ft_prepare_headmodel(cfg_ft, mesh_brain);

    cfg_ft = []; cfg_ft.sourcemodel = sources_brain;
    cfg_ft.headmodel = vol; cfg_ft.grad = braingrad; cfg_ft.reducerank = 'no';
    LF = ft_prepare_leadfield(cfg_ft);

    % LCMV beamformer
    cfg_ft = []; cfg_ft.covariance = 'yes';
    tlock = ft_timelockanalysis(cfg_ft, braindat);

    cfg_ft = []; cfg_ft.method = 'lcmv'; cfg_ft.headmodel = vol;
    cfg_ft.sourcemodel.leadfield = LF; cfg_ft.grid = sources_brain;
    cfg_ft.unit = LF.unit; cfg_ft.lcmv.keepfilter = 'yes';
    source_idx = ft_sourceanalysis(cfg_ft, tlock);

    % Map M1 ROI positions to source grid
    idx_roi = dsearchn(source_idx.pos, intersection_pos);
    idx_roi = unique(idx_roi);
    idx_roi = idx_roi(source_idx.inside(idx_roi));

    roi_center = mean(source_idx.pos(idx_roi,:), 1);
    d = sqrt(sum((intersection_pos - roi_center).^2, 2));
    R = max(d);

    % diagnostic plot (first subject only)
    if ss == 1
        d_all = sqrt(sum((source_idx.pos - roi_center).^2, 2));
        roi_idx_radius = find(d_all <= R & source_idx.inside);
        figure('Color','w');
        ft_plot_mesh(mesh_brain,'facecolor',[0.8 0.8 0.8],'facealpha',0.7,'edgecolor','none');
        hold on;
        plot3(source_idx.pos(roi_idx_radius,1), ...
              source_idx.pos(roi_idx_radius,2), ...
              source_idx.pos(roi_idx_radius,3), ...
              'ro','MarkerSize',6,'MarkerFaceColor','r');
        title(sprintf('Brain VE ROI — %s', sub),'Interpreter','none')
        view(176,-10); camlight;
    end

    % virtual channel via SVD
    cfg_ft = []; cfg_ft.pos = roi_center; cfg_ft.radius = R;
    cfg_ft.method = 'svd'; cfg_ft.numcomponent = 1;
    VE = ft_virtualchannel(cfg_ft, braindat, source_idx);

    savename = sprintf('sub%s_VE_brain_forspectra', sub);
    save(fullfile(save_dir, savename), 'VE')
    fprintf('    Saved: %s\n', savename)
end

end  % run_brain_VE


% -------------------------------------------------------------------------
function run_coherence_spectra(cfg)
% NeuroSpec + FieldTrip pairwise coherence for Brain, Spine, EMG.

subs     = cfg.subs_coh;
save_dir = cfg.save_dir;
data_root = cfg.data_root;
fband    = cfg.fband;
seg_pwr  = cfg.seg_pwr;
lat_min  = cfg.lat_min;
lat_max  = cfg.lat_max;
ve_suffix  = cfg.spine_ve_suffix;
saveFigs   = cfg.saveFigs;
fig_dir    = fullfile(cfg.save_dir, 'figures');
if saveFigs && ~exist(fig_dir,'dir'), mkdir(fig_dir); end

nSubs = numel(subs);

% pre-allocate group arrays
Pstat_brain = nan(nSubs,1);
Prest_brain = nan(nSubs,1);
Pstat_spine = nan(nSubs,1);
Prest_spine = nan(nSubs,1);

results = struct();

for ss = 1:nSubs
    sub = subs{ss};
    fprintf('  [Step B] Subject %s (%d/%d)\n', sub, ss, nSubs);

    datwithEMGmerged = get_datafile(data_root, sub);
    D     = spm_eeg_load(datwithEMGmerged);
    ftdat = spm2fieldtrip(D);

    % load brain VE
    bVE_file = fullfile(save_dir, sprintf('sub%s_VE_brain_forspectra.mat', sub));
    bVE = load(bVE_file); bVE = bVE.VE;

    % load spine VE
    sVE_file = fullfile(save_dir, ...
        sprintf('VE_spine_sub%s_forspectra_bemv2%s.mat', sub, ve_suffix));
    sVE = load(sVE_file); sVE = sVE.VE;

    % EMG (rectified)
    cfg_ft = []; cfg_ft.channel = 'EXG1';
    EMG = ft_selectdata(cfg_ft, ftdat);
    cfg_ft = []; cfg_ft.rectify = 'yes';
    EMG = ft_preprocessing(cfg_ft, EMG);

    % separate conditions
    statidx = find(ftdat.trialinfo==1);
    restidx = find(ftdat.trialinfo==2);
    [nTrials,~] = min([length(statidx) length(restidx)]);

    cfg_ft = []; cfg_ft.trials = statidx(1:nTrials);
    statB   = ft_selectdata(cfg_ft, bVE);
    statS   = ft_selectdata(cfg_ft, sVE);
    statEMG = ft_selectdata(cfg_ft, EMG);

    cfg_ft.trials = restidx(1:nTrials);
    restB   = ft_selectdata(cfg_ft, bVE);
    restS   = ft_selectdata(cfg_ft, sVE);
    restEMG = ft_selectdata(cfg_ft, EMG);

    %% ---- FieldTrip coherence spectra (static condition) ----
    statB.label{1}   = 'brain';
    statS.label{1}   = 'spine';
    statEMG.label{1} = 'EMG';
    alldat = ft_appenddata([], statB, statS, statEMG);

    cfg_ft = []; cfg_ft.output = 'fourier'; cfg_ft.method = 'mtmfft';
    cfg_ft.foilim = [2 75]; cfg_ft.tapsmofrq = 2; cfg_ft.keeptrials = 'yes';
    freq = ft_freqanalysis(cfg_ft, alldat);

    cfg_ft = []; cfg_ft.method = 'coh';
    coh = ft_connectivityanalysis(cfg_ft, freq);

    %% ---- Spectral power: static vs rest ----
    statB2 = statB; statB2.label{1} = 'brain';
    statS2 = statS; statS2.label{1} = 'spine';
    restB2 = restB; restB2.label{1} = 'brain';
    restS2 = restS; restS2.label{1} = 'spine';

    statBS = ft_appenddata([], statB2, statS2);
    restBS = ft_appenddata([], restB2, restS2);

    cfgp = []; cfgp.output = 'pow'; cfgp.method = 'mtmfft';
    cfgp.foilim = [2 75]; cfgp.tapsmofrq = 2; cfgp.keeptrials = 'yes';
    freq_stat = ft_freqanalysis(cfgp, statBS);
    freq_rest = ft_freqanalysis(cfgp, restBS);

    % log power spectrum plot (subject 1 only)
    if ss == 1
        iSp = find(strcmp(freq_stat.label,'spine'));
        Psp_stat = squeeze(mean(freq_stat.powspctrm(:,iSp,:),1));
        Psp_rest = squeeze(mean(freq_rest.powspctrm(:,iSp,:),1));

        hfig_pow1 = figure('Color','w'); hold on;
        plot(freq_stat.freq, log10(Psp_stat), 'LineWidth',2);
        plot(freq_rest.freq, log10(Psp_rest), 'LineWidth',2);
        xlim([2 75]); xlabel('Frequency (Hz)'); ylabel('log_{10} Power');
        title(sprintf('Participant 1 (%s): Spinal VE log power', sub),'Interpreter','none');
        legend({'Contraction','Rest'},'Location','best'); box off;
        if saveFigs, savefig(hfig_pow1, fullfile(fig_dir, sprintf('P1_%s_spinalVE_logpower.fig', sub))); end

        % band power t-test for subject 1
        cfgsel = []; cfgsel.frequency = fband; cfgsel.avgoverfreq = 'yes';
        bp_stat1 = ft_selectdata(cfgsel, freq_stat);
        bp_rest1 = ft_selectdata(cfgsel, freq_rest);
        x_stat1  = log(squeeze(bp_stat1.powspctrm));
        x_rest1  = log(squeeze(bp_rest1.powspctrm));
        iB1 = find(strcmp(bp_stat1.label,'brain'));
        iS1 = find(strcmp(bp_stat1.label,'spine'));

        [~,pB,~,sB] = ttest2(x_stat1(:,iB1), x_rest1(:,iB1), 'Vartype','unequal');
        [~,pS,~,sS] = ttest2(x_stat1(:,iS1), x_rest1(:,iS1), 'Vartype','unequal');
        sdB = std([x_stat1(:,iB1); x_rest1(:,iB1)]);
        sdS = std([x_stat1(:,iS1); x_rest1(:,iS1)]);
        fprintf('\n  Participant 1 band power (10-35 Hz, log) Static vs Rest:\n');
        fprintf('    Brain: t=%.3f, p=%.4g, d=%.3f\n', sB.tstat, pB, ...
            (mean(x_stat1(:,iB1))-mean(x_rest1(:,iB1)))/sdB);
        fprintf('    Spine: t=%.3f, p=%.4g, d=%.3f\n', sS.tstat, pS, ...
            (mean(x_stat1(:,iS1))-mean(x_rest1(:,iS1)))/sdS);
    end

    % average band power per subject
    cfgsel = []; cfgsel.frequency = fband; cfgsel.avgoverfreq = 'yes';
    bp_stat = ft_selectdata(cfgsel, freq_stat);
    bp_rest = ft_selectdata(cfgsel, freq_rest);
    m_stat  = squeeze(mean(bp_stat.powspctrm, 1));
    m_rest  = squeeze(mean(bp_rest.powspctrm, 1));
    iB = find(strcmp(bp_stat.label,'brain'));
    iS = find(strcmp(bp_stat.label,'spine'));
    Pstat_brain(ss) = m_stat(iB);
    Prest_brain(ss) = m_rest(iB);
    Pstat_spine(ss) = m_stat(iS);
    Prest_spine(ss) = m_rest(iS);

    %% ---- NeuroSpec (static condition only) ----
    statBcont   = [statB.trial{:}];
    statScont   = [statS.trial{:}];
    statEMGcont = abs([statEMG.trial{:}]);   % rectify

    samp_rate = ftdat.hdr.Fs;
    opt_str   = 'M3';

    [f1,t1,~] = sp2a2_R2_mt(statBcont', statEMGcont', samp_rate, seg_pwr);
    [f2,t2,~] = sp2a2_R2_mt(statBcont', statScont',   samp_rate, seg_pwr, opt_str);
    [f3,t3,~] = sp2a2_R2_mt(statScont', statEMGcont', samp_rate, seg_pwr, opt_str);

    % partial coherence brain-EMG conditioned on spine
    [fp,~,~] = sp2a2_R2_pc1(statBcont', statEMGcont', statScont', samp_rate, seg_pwr);

    % %% ---- Individual coherence plot (3 pairs) ----
    % freq_band = fband;
    % plot_colored_band = @(f, y, band, clr) ...
    %     plot(f(f>=band(1) & f<=band(2)), y(f>=band(1) & f<=band(2)), ...
    %          'Color', clr, 'LineWidth', 2);
% 
    % hfig_coh = figure('Color','w','Position',[100 100 600 800]);
% 
    % subplot(3,1,1); hold on;
    % plot(f1(:,1), f1(:,12), 'Color',[0.8 0.8 0.8], 'LineWidth',2);
    % plot(f1(:,1), f1(:,11), 'Color',[0.8 0.8 0.8], 'LineWidth',2);
    % hRev = plot_colored_band(f1(:,1), f1(:,12), freq_band, [0.2 0.4 0.8]);
    % hFwd = plot_colored_band(f1(:,1), f1(:,11), freq_band, [0.8 0.2 0.2]);
    % xlim([0 40]); xlabel('Frequency (Hz)','FontSize',14);
    % ylabel('Coherence','FontSize',14); title('Brain \leftrightarrow EMG','FontSize',16);
    % legend([hRev hFwd],{'Reverse','Forward'},'Location','northeast'); box off;
% 
    % subplot(3,1,2); hold on;
    % plot(f2(:,1), f2(:,12), 'Color',[0.8 0.8 0.8], 'LineWidth',2);
    % plot(f2(:,1), f2(:,11), 'Color',[0.8 0.8 0.8], 'LineWidth',2);
    % hRev = plot_colored_band(f2(:,1), f2(:,12), freq_band, [0.2 0.4 0.8]);
    % hFwd = plot_colored_band(f2(:,1), f2(:,11), freq_band, [0.8 0.2 0.2]);
    % xlim([0 40]); xlabel('Frequency (Hz)','FontSize',14);
    % ylabel('Coherence','FontSize',14); title('Brain \leftrightarrow Spine','FontSize',16);
    % legend([hRev hFwd],{'Reverse','Forward'},'Location','northeast'); box off;
% 
    % subplot(3,1,3); hold on;
    % plot(f3(:,1), f3(:,12), 'Color',[0.8 0.8 0.8], 'LineWidth',2);
    % plot(f3(:,1), f3(:,11), 'Color',[0.8 0.8 0.8], 'LineWidth',2);
    % hRev = plot_colored_band(f3(:,1), f3(:,12), freq_band, [0.2 0.4 0.8]);
    % hFwd = plot_colored_band(f3(:,1), f3(:,11), freq_band, [0.8 0.2 0.2]);
    % xlim([0 40]); xlabel('Frequency (Hz)','FontSize',14);
    % ylabel('Coherence','FontSize',14); title('Spine \leftrightarrow EMG','FontSize',16);
    % legend([hRev hFwd],{'Reverse','Forward'},'Location','northeast'); box off;
% 
    % sgtitle(sprintf('Subject %s — Directed coherence (static)', sub), 'Interpreter','none');
    % if saveFigs, savefig(hfig_coh, fullfile(fig_dir, sprintf('sub%s_directed_coherence.fig', sub))); end
% 
    % % partial coherence overlay (brain-EMG vs partial)
    % hfig_part = figure('Color','w'); hold on;
    % plot(fp(:,1), fp(:,4), 'LineWidth',2);
    % plot(f1(:,1), f1(:,4),  'LineWidth',2);
    % legend({'Partial (conditioned on spine)','Full'},'Location','best');
    % xlim([0 45]); xlabel('Frequency (Hz)'); ylabel('Coherence');
    % title(sprintf('%s: Brain-EMG coherence (full vs partial)', sub),'Interpreter','none');
    % box off;
    % if saveFigs, savefig(hfig_part, fullfile(fig_dir, sprintf('sub%s_brainEMG_partial_coherence.fig', sub))); end

    %% ---- Directed coherence areas ----
    freq_band = fband;
    compute_area = @(fmat, fwd_col, rev_col, band) deal( ...
        trapz(fmat(fmat(:,1)>=band(1) & fmat(:,1)<=band(2),1), ...
              fmat(fmat(:,1)>=band(1) & fmat(:,1)<=band(2),fwd_col)), ...
        trapz(fmat(fmat(:,1)>=band(1) & fmat(:,1)<=band(2),1), ...
              fmat(fmat(:,1)>=band(1) & fmat(:,1)<=band(2),rev_col)) );

    [brainEMG_fwd,   brainEMG_rev]   = compute_area(f1, 11, 12, freq_band);
    [brainSpine_fwd, brainSpine_rev] = compute_area(f2, 11, 12, freq_band);
    [spineEMG_fwd,   spineEMG_rev]   = compute_area(f3, 11, 12, freq_band);

    results(ss).sub = sub;
    results(ss).brainEMG.forward_area   = brainEMG_fwd;
    results(ss).brainEMG.reverse_area   = brainEMG_rev;
    results(ss).brainEMG.ratio          = brainEMG_fwd / (brainEMG_fwd + brainEMG_rev);
    results(ss).brainSpine.forward_area = brainSpine_fwd;
    results(ss).brainSpine.reverse_area = brainSpine_rev;
    results(ss).brainSpine.ratio        = brainSpine_fwd / (brainSpine_fwd + brainSpine_rev);
    results(ss).spineEMG.forward_area   = spineEMG_fwd;
    results(ss).spineEMG.reverse_area   = spineEMG_rev;
    results(ss).spineEMG.ratio          = spineEMG_fwd / (spineEMG_fwd + spineEMG_rev);

    %% ---- Peak cross-correlation latencies ----
    [results(ss).brainEMG.peak_rho,   results(ss).brainEMG.peak_latency]   = ...
        select_peak(t1, lat_min, lat_max);
    [results(ss).brainSpine.peak_rho, results(ss).brainSpine.peak_latency] = ...
        select_peak(t2, lat_min, lat_max);
    [results(ss).spineEMG.peak_rho,   results(ss).spineEMG.peak_latency]   = ...
        select_peak(t3, lat_min, lat_max);

    % store top-3 peaks (abs latency) for distribution plot
    results(ss).brainEMG.top_lats   = get_top_peaks(t1, lat_min, lat_max, 3);
    results(ss).brainSpine.top_lats = get_top_peaks(t2, lat_min, lat_max, 3);
    results(ss).spineEMG.top_lats   = get_top_peaks(t3, lat_min, lat_max, 3);

end  % subject loop

%% ---- Group: paired t-tests band power ----
fprintf('\n=== Group: paired t-tests 10-35 Hz band power (Static vs Rest) ===\n');

okB = isfinite(Pstat_brain) & isfinite(Prest_brain);
[~,pb,~,sb] = ttest(log(Pstat_brain(okB)), log(Prest_brain(okB)));
fprintf('  Brain: n=%d, t(%d)=%.3f, p=%.4g\n', sum(okB), sb.df, sb.tstat, pb);

okS = isfinite(Pstat_spine) & isfinite(Prest_spine);
[~,ps,~,ss_] = ttest(log(Pstat_spine(okS)), log(Prest_spine(okS)));
fprintf('  Spine: n=%d, t(%d)=%.3f, p=%.4g\n', sum(okS), ss_.df, ss_.tstat, ps);

%% ---- Group: directed coherence bar plot ----
forward_areas = zeros(nSubs,3);
reverse_areas = zeros(nSubs,3);
for ss = 1:nSubs
    forward_areas(ss,1) = results(ss).brainEMG.forward_area;
    forward_areas(ss,2) = results(ss).brainSpine.forward_area;
    forward_areas(ss,3) = results(ss).spineEMG.forward_area;
    reverse_areas(ss,1) = results(ss).brainEMG.reverse_area;
    reverse_areas(ss,2) = results(ss).brainSpine.reverse_area;
    reverse_areas(ss,3) = results(ss).spineEMG.reverse_area;
end

mean_fwd = mean(forward_areas,1);
mean_rev = mean(reverse_areas,1);
sem_fwd  = std(forward_areas,0,1)/sqrt(nSubs);
sem_rev  = std(reverse_areas,0,1)/sqrt(nSubs);

data_bar = [mean_fwd; mean_rev]';
sems_bar = [sem_fwd;  sem_rev]';
x_labels = {'Brain\rightarrowEMG','Brain\rightarrowSpine','Spine\rightarrowEMG'};

hfig_bar = figure('Color','w'); hold on;
b = bar(data_bar,'grouped');
b(1).FaceColor = [0.3 0.6 0.9]; b(1).FaceAlpha = 0.5;
b(2).FaceColor = [0.9 0.3 0.3]; b(2).FaceAlpha = 0.5;

[ngroups, nbars] = size(data_bar);
groupwidth = min(0.8, nbars/(nbars+1.5));
for i = 1:nbars
    xpos = (1:ngroups) - groupwidth/2 + (2*i-1)*groupwidth/(2*nbars);
    errorbar(xpos, data_bar(:,i), sems_bar(:,i), 'k.', 'LineWidth',1.5);
end

% individual points + participant 1 circle
p1 = 1;
jitter = 0.015;
for i = 1:nbars
    xpos = (1:ngroups) - groupwidth/2 + (2*i-1)*groupwidth/(2*nbars);
    for j = 1:ngroups
        y  = forward_areas(:,j);
        if i == 2, y = reverse_areas(:,j); end
        xj = xpos(j) + jitter*randn(nSubs,1);
        scatter(xj, y, 45, 'MarkerFaceColor',[0.35 0.35 0.35], ...
            'MarkerEdgeColor','none','MarkerFaceAlpha',0.55);
        scatter(xj(p1), y(p1), 180, 'o', 'MarkerEdgeColor','k', ...
            'LineWidth',2,'MarkerFaceColor','none');
    end
end

set(gca,'XTick',1:ngroups,'XTickLabel',x_labels,'FontSize',14);
ylabel('Coherence area (10-35 Hz)');
ylim([0 max(data_bar(:))*1.3]);
legend({'Forward','Reverse'},'Location','northwest'); box off;
title('Group directed coherence','Interpreter','none')
if saveFigs, savefig(hfig_bar, fullfile(fig_dir, 'group_directed_coherence.fig')); end

%% ---- Group: peak latency plot ----
brainEMG_lat   = abs(arrayfun(@(s) s.brainEMG.peak_latency,   results));
brainSpine_lat = abs(arrayfun(@(s) s.brainSpine.peak_latency, results));
spineEMG_lat   = abs(arrayfun(@(s) s.spineEMG.peak_latency,   results));

lightGrey  = [0.75 0.75 0.75];
darkGrey   = [0.25 0.25 0.25];
medianGrey = [0.10 0.10 0.10];
x_lat      = [1 2 3];
p1         = 1;

hfig_lat = figure('Color','w'); hold on;
hP1  = [];
hMed = [];

% individual lines
for ss = 1:nSubs
    y = [brainEMG_lat(ss), brainSpine_lat(ss), spineEMG_lat(ss)];
    if any(isnan(y)), continue; end
    if ss == p1
        hP1 = plot(x_lat, y, '-o', ...
            'LineWidth',3.5,'MarkerSize',9, ...
            'Color',darkGrey,'MarkerFaceColor',darkGrey,'MarkerEdgeColor','w');
    else
        plot(x_lat, y, '-','Color',lightGrey,'LineWidth',0.8);
    end
end

% scatter points
scatter(ones(nSubs,1),   brainEMG_lat',   60, lightGrey, 'filled','MarkerFaceAlpha',0.5,'MarkerEdgeColor','none');
scatter(2*ones(nSubs,1), brainSpine_lat', 60, lightGrey, 'filled','MarkerFaceAlpha',0.5,'MarkerEdgeColor','none');
scatter(3*ones(nSubs,1), spineEMG_lat',   60, lightGrey, 'filled','MarkerFaceAlpha',0.5,'MarkerEdgeColor','none');

% median + MAD
med_vals = [median(brainEMG_lat,'omitnan'), ...
            median(brainSpine_lat,'omitnan'), ...
            median(spineEMG_lat,'omitnan')];
mad_vals = [mad(brainEMG_lat,1), mad(brainSpine_lat,1), mad(spineEMG_lat,1)];

hMed = plot(x_lat, med_vals, '--s', ...
    'LineWidth',2.5,'MarkerSize',9, ...
    'Color',medianGrey,'MarkerFaceColor','w','MarkerEdgeColor',medianGrey);
errorbar(x_lat, med_vals, mad_vals, ...
    'Color',medianGrey,'LineStyle','none','LineWidth',1.5,'CapSize',10);

% legend
leg_handles = [];
leg_labels  = {};
if ~isempty(hP1)
    leg_handles = [leg_handles, hP1];
    leg_labels  = [leg_labels, {'Participant 1'}];
end
leg_handles = [leg_handles, hMed];
leg_labels  = [leg_labels,  {'Median +/- MAD'}];
legend(leg_handles, leg_labels, 'Location','northwest','Box','off');

xlim([0.5 3.5]); xticks(x_lat);
xticklabels({'Brain<->EMG','Brain<->Spine','Spine<->EMG'});
ylabel('Peak latency (ms)'); set(gca,'FontSize',14); grid on; box on;
title('Cross-correlation peak latencies','Interpreter','none')
if saveFigs, savefig(hfig_lat, fullfile(fig_dir, 'group_peak_latencies.fig')); end


%% ---- Print latency summary ----
fprintf('\n=== Peak cross-correlation latencies (ms) ===\n');
fprintf('  Pair            P1      Median   MAD\n');
fprintf('  Brain<->EMG   %5.1f   %5.1f    %5.1f\n', ...
    brainEMG_lat(p1), med_vals(1), mad_vals(1));
fprintf('  Brain<->Spine %5.1f   %5.1f    %5.1f\n', ...
    brainSpine_lat(p1), med_vals(2), mad_vals(2));
fprintf('  Spine<->EMG   %5.1f   %5.1f    %5.1f\n', ...
    spineEMG_lat(p1), med_vals(3), mad_vals(3));
fprintf('=============================================\n');

end  % run_coherence_spectra


%% =========================================================================
%  SHARED UTILITY FUNCTIONS
%% =========================================================================

function fname = get_datafile(data_root, sub)
    run = '001'; if strcmp(sub,'OP00224'), run = '002'; end
    fname = fullfile(data_root, ['sub-' sub], 'ses-001', 'meg', ...
        sprintf('pmergedoe1000mspddfflo45hi45hfcstatic_%s_array1.mat', run));
end

function lats = get_top_peaks(t, lat_min, lat_max, nPeaks)
% Return abs latencies of top N peaks by |rho| in the given window.
    valid_idx = t(:,1) >= lat_min & t(:,1) <= lat_max;
    rho_vals  = t(valid_idx, 3);
    lag_vals  = t(valid_idx, 1);

    [~, locs] = findpeaks(abs(rho_vals));
    if isempty(locs), locs = (1:numel(rho_vals))'; end
    locs = locs(abs(lag_vals(locs)) > 2);   % exclude ±2 ms
    if isempty(locs)
        [~, locs] = max(abs(rho_vals));
    end

    [~, ord] = sort(abs(rho_vals(locs)), 'descend');
    locs_sorted = locs(ord);
    nOut = min(nPeaks, numel(locs_sorted));
    lats = lag_vals(locs_sorted(1:nOut));
end

function braingrad = subset_grad(grad_mm, idx)
    braingrad          = grad_mm;
    braingrad.chanori  = grad_mm.chanori(idx,:);
    braingrad.chanpos  = grad_mm.chanpos(idx,:);
    braingrad.chantype = grad_mm.chantype(1:length(idx));
    braingrad.chanunit = grad_mm.chanunit(1:length(idx));
    braingrad.coilori  = grad_mm.coilori(idx,:);
    braingrad.coilpos  = grad_mm.coilpos(idx,:);
    braingrad.label    = grad_mm.label(idx);
    braingrad.tra      = grad_mm.tra(idx,idx);
end
