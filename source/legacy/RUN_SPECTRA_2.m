%% RUN_SPECTRA.m
% Downstream spectra and coherence analysis for brain-spine pipeline.
%
% REQUIRES outputs from RUN_PIPELINE.m:
%   - VE_spine_sub*_bemv2*.mat    (Step 3) - spinal virtual electrodes
%   - sub*_VE_brain_forspectra.mat          - brain virtual electrodes (legacy)
%
% STEPS:
%   Step B - Pairwise coherence spectra (NeuroSpec): Brain-EMG, Brain-Spine, Spine-EMG
%             + directed coherence (forward/reverse areas)
%             + FieldTrip coherence spectra
%             + spectral power (10-35 Hz) static vs rest
%             + group coherence spectra (contraction vs rest)
%             + peak frequency extraction per pair (low/high beta)
%             + partial coherence brain-EMG conditioned on spine
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
cfg.data_root = 'C:\spinecoh_data';
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
cfg.doSmooth             = 1;
cfg.spine_smooth_fwhm_mm = 20;   % must match RUN_PIPELINE

% --- Analysis parameters --------------------------------------------------
cfg.fband    = [10 35];   % beta band for power/coherence summary
cfg.seg_pwr  = 11;        % NeuroSpec segment power (2^11 = 2048 samples)
cfg.lat_min  = -50;       % cross-correlation latency window (ms)
cfg.lat_max  =  50;

% --- Figure saving -------------------------------------------------------
saveFigs = 1;   % 1 = save all figures as .fig to save_dir/figures/

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

fprintf('\n=== SPECTRA CONFIG ===\n')
fprintf('  Brain subs: %d  |  Spine subs: %d  |  Coh subs: %d\n', ...
    numel(cfg.subs_brain), numel(cfg.subs_spine), numel(cfg.subs_coh))
fprintf('======================\n\n')

%% =========================================================================
%  STEP B - Coherence spectra (NeuroSpec + FieldTrip)
%% =========================================================================

fprintf('\n>>> STEP B: Coherence spectra\n')
run_coherence_spectra(cfg)
fprintf('>>> STEP B complete.\n\n')

fprintf('\n=== SPECTRA FINISHED ===\n\n')

% -------------------------------------------------------------------------
function run_coherence_spectra(cfg)
% NeuroSpec + FieldTrip pairwise coherence for Brain, Spine, EMG.

subs      = cfg.subs_coh;
save_dir  = cfg.save_dir;
data_root = cfg.data_root;
fband     = cfg.fband;
seg_pwr   = cfg.seg_pwr;
lat_min   = cfg.lat_min;
lat_max   = cfg.lat_max;
ve_suffix = cfg.spine_ve_suffix;
saveFigs  = cfg.saveFigs;
fig_dir   = fullfile(cfg.save_dir, 'figures');
if saveFigs && ~exist(fig_dir,'dir'), mkdir(fig_dir); end

nSubs  = numel(subs);
nFreqs = [];
all_dir_brainEMG_stat   = [];
all_dir_brainSpine_stat = [];
all_dir_spineEMG_stat   = [];

peak_freq_brainEMG   = nan(nSubs,1);
peak_freq_brainSpine = nan(nSubs,1);
peak_freq_spineEMG   = nan(nSubs,1);
peak_coh_brainEMG    = nan(nSubs,1);
peak_coh_brainSpine  = nan(nSubs,1);
peak_coh_spineEMG    = nan(nSubs,1);
peak_coh_brainEMG_nt   = nan(nSubs,1);
peak_coh_brainSpine_nt = nan(nSubs,1);
peak_coh_spineEMG_nt   = nan(nSubs,1);

freq_axis = [];

Pstat_brain = nan(nSubs,1);
Prest_brain = nan(nSubs,1);
Pstat_spine = nan(nSubs,1);
Prest_spine = nan(nSubs,1);
Pstat_brain_fooof = nan(nSubs,1);
Prest_brain_fooof = nan(nSubs,1);
Pstat_spine_fooof = nan(nSubs,1);
Prest_spine_fooof = nan(nSubs,1);

results = struct();
psi_brainEMG   = nan(nSubs,1);
psi_brainSpine = nan(nSubs,1);
psi_spineEMG   = nan(nSubs,1);

col_BE = [0.2 0.4 0.8];
col_BS = [0.8 0.2 0.4];
col_SE = [0.2 0.7 0.3];

peak_win_hz = 3;
M_segments  = nan;

for ss = 1:nSubs
    sub = subs{ss};
    fprintf('  [Step B] Subject %s (%d/%d)\n', sub, ss, nSubs);

    datwithEMGmerged = get_datafile(data_root, sub);
    D     = spm_eeg_load(datwithEMGmerged);
    ftdat = spm2fieldtrip(D);

    % load brain VE (legacy)
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

    %% ---- FieldTrip coherence spectra (contraction condition) ----
    statB.label{1}   = 'brain';
    statS.label{1}   = 'spine';
    statEMG.label{1} = 'EMG';
    alldat = ft_appenddata([], statB, statS, statEMG);

    cfg_ft = []; cfg_ft.output = 'fourier'; cfg_ft.method = 'mtmfft';
    cfg_ft.foilim = [2 75]; cfg_ft.tapsmofrq = 2; cfg_ft.keeptrials = 'yes';
    freq = ft_freqanalysis(cfg_ft, alldat);

    cfg_ft = []; cfg_ft.method = 'coh';
    coh = ft_connectivityanalysis(cfg_ft, freq);

    %% ---- Spectral power: contraction vs rest ----
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

    % FOOOF periodic component
    try
        cfgf = []; cfgf.output = 'fooof_peaks'; cfgf.method = 'mtmfft';
        cfgf.foilim = [2 75]; cfgf.tapsmofrq = 2; cfgf.keeptrials = 'no';
        freq_stat_fooof = ft_freqanalysis(cfgf, statBS);
        freq_rest_fooof = ft_freqanalysis(cfgf, restBS);

        cfgsel_f = []; cfgsel_f.frequency = fband; cfgsel_f.avgoverfreq = 'yes';
        bp_stat_f = ft_selectdata(cfgsel_f, freq_stat_fooof);
        bp_rest_f = ft_selectdata(cfgsel_f, freq_rest_fooof);
        m_stat_f = squeeze(bp_stat_f.powspctrm);
        m_rest_f = squeeze(bp_rest_f.powspctrm);
        iBf = find(strcmp(bp_stat_f.label,'brain'));
        iSf = find(strcmp(bp_stat_f.label,'spine'));
        Pstat_brain_fooof(ss) = m_stat_f(iBf);
        Prest_brain_fooof(ss) = m_rest_f(iBf);
        Pstat_spine_fooof(ss) = m_stat_f(iSf);
        Prest_spine_fooof(ss) = m_rest_f(iSf);

        if ss == 1
            fooof_stat_P1 = freq_stat_fooof;
            fooof_rest_P1 = freq_rest_fooof;
        end
    catch ME
        warning('FOOOF failed for subject %s: %s', sub, ME.message);
    end

    %% ---- NeuroSpec (contraction) ----
    statBcont   = [statB.trial{:}];
    statScont   = [statS.trial{:}];
    statEMGcont = abs([statEMG.trial{:}]);

    samp_rate = ftdat.hdr.Fs;
    opt_str   = 'M4';

    [f1_nt, t1_nt, ~] = sp2a2_R2_mt(statBcont', statEMGcont', samp_rate, seg_pwr);
    [f2_nt, t2_nt, ~] = sp2a2_R2_mt(statBcont', statScont',   samp_rate, seg_pwr);
    [f3_nt, t3_nt, ~] = sp2a2_R2_mt(statScont', statEMGcont', samp_rate, seg_pwr);
    [f1_mt, t1_mt, ~] = sp2a2_R2_mt(statBcont', statEMGcont', samp_rate, seg_pwr, opt_str);
    [f2_mt, t2_mt, ~] = sp2a2_R2_mt(statBcont', statScont',   samp_rate, seg_pwr, opt_str);
    [f3_mt, t3_mt, ~] = sp2a2_R2_mt(statScont', statEMGcont', samp_rate, seg_pwr, opt_str);

    % Partial coherence - Brain-EMG conditioned on spine (MT)
    mt_NW      = 4;
    mt_seg_len = 2^seg_pwr;
    [fp_mt_BE_spine, freq_mt, fc_mt_BE] = mt_partial_coherence(statBcont', statEMGcont', statScont', samp_rate, mt_seg_len, mt_NW);
    [~, ~,        fc_mt_SE] = mt_partial_coherence(statScont', statEMGcont', statBcont', samp_rate, mt_seg_len, mt_NW);
    [~,   ~,      fc_mt_BS] = mt_partial_coherence(statBcont', statScont',   statEMGcont', samp_rate, mt_seg_len, mt_NW);

    ns_len = size(f1_mt, 1);
    fp_mt_BE_spine = fp_mt_BE_spine(1:ns_len);
    fc_mt_BE       = fc_mt_BE(1:ns_len);
    fc_mt_SE       = fc_mt_SE(1:ns_len);
    fc_mt_BS       = fc_mt_BS(1:ns_len);
    freq_mt        = freq_mt(1:ns_len);

    if isnan(M_segments)
        seg_len    = 2^seg_pwr;
        M_segments = floor(length(statBcont) / seg_len);
        fprintf('  NeuroSpec: %d segments of length %d (fs=%d Hz)\n', ...
            M_segments, seg_len, samp_rate);
    end

    if isempty(nFreqs)
        nFreqs = size(f1_mt, 1);
        all_dir_brainEMG_stat   = nan(nSubs, nFreqs);
        all_dir_brainSpine_stat = nan(nSubs, nFreqs);
        all_dir_spineEMG_stat   = nan(nSubs, nFreqs);
        fprintf('  Array sized to %d frequency bins from first subject output.\n', nFreqs);
    end

    freq_axis = f1_mt(:,1)';
    all_dir_brainEMG_stat(ss,:)   = (f1_mt(:,11) + f1_mt(:,12))';
    all_dir_brainSpine_stat(ss,:) = (f2_mt(:,11) + f2_mt(:,12))';
    all_dir_spineEMG_stat(ss,:)   = (f3_mt(:,11) + f3_mt(:,12))';

    band_mask    = freq_axis >= fband(1) & freq_axis <= fband(2);
    freqs_inband = freq_axis(band_mask);

    [pk1, i1] = max(f1_mt(band_mask, 4));
    peak_freq_brainEMG(ss)  = freqs_inband(i1);
    peak_coh_brainEMG(ss)   = pk1;

    [pk2, i2] = max(f2_mt(band_mask, 4));
    peak_freq_brainSpine(ss) = freqs_inband(i2);
    peak_coh_brainSpine(ss)  = pk2;

    [pk3, i3] = max(f3_mt(band_mask, 4));
    peak_freq_spineEMG(ss)  = freqs_inband(i3);
    peak_coh_spineEMG(ss)   = pk3;

    % Sub-band peak frequencies (low/high beta) - filled after beta_split_hz known
    % Stored using directed sum (col 11+12) for consistency with plots
    dir_BE = f1_mt(:,11) + f1_mt(:,12);
    dir_BS = f2_mt(:,11) + f2_mt(:,12);
    dir_SE = f3_mt(:,11) + f3_mt(:,12);
    results(ss).dir_BE = dir_BE';
    results(ss).dir_BS = dir_BS';
    results(ss).dir_SE = dir_SE';

    fax_nt_pk = f1_nt(:,1)';
    band_nt   = fax_nt_pk >= fband(1) & fax_nt_pk <= fband(2);
    peak_coh_brainEMG_nt(ss)   = max(f1_nt(band_nt, 4));
    peak_coh_brainSpine_nt(ss) = max(f2_mt(band_nt, 4));
    peak_coh_spineEMG_nt(ss)   = max(f3_mt(band_nt, 4));

    % No-taper peak FREQUENCY (for MT vs no-taper comparison)
    freqs_nt_inband = fax_nt_pk(band_nt);
    [~, i1_nt] = max(f1_nt(band_nt, 4));
    [~, i2_nt] = max(f2_mt(band_nt, 4));
    [~, i3_nt] = max(f3_mt(band_nt, 4));
    results(ss).peak_freq_BE_nt  = freqs_nt_inband(i1_nt);
    results(ss).peak_freq_BS_nt  = freqs_nt_inband(i2_nt);
    results(ss).peak_freq_SE_nt  = freqs_nt_inband(i3_nt);

    %% ---- Individual subject: coherence spectra plot (OP00212 only) ----
    if strcmp(sub, 'OP00212')

        hfig_coh = figure('Color','w','Position',[100 100 500 750]);

        % subplot 1: Brain-EMG
        subplot(3,1,1); hold on;
        dir1 = (f1_mt(:,11) + f1_mt(:,12)) * 1e3;
        yl1 = [0 max(dir1) * 1.3];
        plot(f1_mt(:,1), dir1, 'Color',col_BE, 'LineWidth',2);
        xline(fband(1),'--k','Alpha',0.3); xline(fband(2),'--k','Alpha',0.3);
        xlim([2 45]); ylim(yl1);
        xlabel('Frequency (Hz)','FontSize',14);
        ylabel('Coherence (x1e-3)','FontSize',14);
        title('Brain - EMG','FontSize',14);
        set(gca,'FontSize',14); box off;

        % subplot 2: Brain-Spine
        subplot(3,1,2); hold on;
        dir2 = (f2_mt(:,11) + f2_mt(:,12)) * 1e3;
        yl2 = [0 max(dir2) * 1.3];
        plot(f2_mt(:,1), dir2, 'Color',col_BS, 'LineWidth',2);
        xline(fband(1),'--k','Alpha',0.3); xline(fband(2),'--k','Alpha',0.3);
        xlim([2 45]); ylim(yl2);
        xlabel('Frequency (Hz)','FontSize',14);
        ylabel('Coherence (x1e-3)','FontSize',14);
        title('Brain - Spine','FontSize',14);
        set(gca,'FontSize',14); box off;

        % subplot 3: Spine-EMG
        subplot(3,1,3); hold on;
        dir3 = (f3_mt(:,11) + f3_mt(:,12)) * 1e3;
        yl3 = [0 max(dir3) * 1.3];
        plot(f3_mt(:,1), dir3, 'Color',col_SE, 'LineWidth',2);
        xline(fband(1),'--k','Alpha',0.3); xline(fband(2),'--k','Alpha',0.3);
        xlim([2 45]); ylim(yl3);
        xlabel('Frequency (Hz)','FontSize',14);
        ylabel('Coherence (x1e-3)','FontSize',14);
        title('Spine - EMG','FontSize',14);
        set(gca,'FontSize',14); box off;

        sgtitle(sprintf('Subject %s - non-zero lag coherence', sub), ...
            'Interpreter','none','FontSize',14);
        if saveFigs
            savefig(hfig_coh, fullfile(fig_dir, ...
                sprintf('sub%s_coherence_spectra.fig', sub)));
        end

        % Separate figure: Brain-EMG full vs partial (conditioned on spine)
        fax_mt_sub = f1_mt(:,1);
        ylp = [0 max([fc_mt_BE; fp_mt_BE_spine])*1e3*1.3];
        hfig_part_p1 = figure('Color','w','Position',[100 100 500 350]);
        hold on;
        hfull_s = plot(fax_mt_sub, fc_mt_BE*1e3,      'Color', col_BE,        'LineWidth', 2);
        hpart_s = plot(fax_mt_sub, fp_mt_BE_spine*1e3, 'Color', [0.6 0.2 0.8], 'LineWidth', 2, 'LineStyle', '--');
        xline(fband(1),'--k','Alpha',0.3); xline(fband(2),'--k','Alpha',0.3);
        xlim([2 45]); ylim(ylp);
        xlabel('Frequency (Hz)','FontSize',14); ylabel('Coherence (x1e-3)','FontSize',14);
        title(sprintf('Subject %s - Brain-EMG | spine (MT)', sub), 'Interpreter','none','FontSize',14);
        legend([hfull_s hpart_s],{'Full','Partial (spine removed)'},'Location','northeast','Box','off','FontSize',12);
        set(gca,'FontSize',14); box off;
        if saveFigs
            savefig(hfig_part_p1, fullfile(fig_dir, ...
                sprintf('sub%s_partial_coherence_BE_spine.fig', sub)));
        end
    end  % OP00212 only

    %% ---- Directed coherence areas (no-taper for all pairs) ----
    freq_band = fband;
    compute_area = @(fmat, fwd_col, rev_col, band) deal( ...
        trapz(fmat(fmat(:,1)>=band(1) & fmat(:,1)<=band(2),1), ...
              fmat(fmat(:,1)>=band(1) & fmat(:,1)<=band(2),fwd_col)), ...
        trapz(fmat(fmat(:,1)>=band(1) & fmat(:,1)<=band(2),1), ...
              fmat(fmat(:,1)>=band(1) & fmat(:,1)<=band(2),rev_col)) );

    [brainEMG_fwd,   brainEMG_rev]   = compute_area(f1_nt, 11, 12, freq_band);
    [brainSpine_fwd, brainSpine_rev] = compute_area(f2_nt, 11, 12, freq_band);
    [spineEMG_fwd,   spineEMG_rev]   = compute_area(f3_nt, 11, 12, freq_band);

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

    results(ss).partial_mt_BE_spine = fp_mt_BE_spine';
    results(ss).full_coh_mt_BE      = fc_mt_BE';
    results(ss).full_coh_mt_SE      = fc_mt_SE';
    results(ss).full_coh_mt_BS      = fc_mt_BS';
    results(ss).freq_axis_mt        = f1_mt(:,1)';
    results(ss).peak_freq_BE_mt     = peak_freq_brainEMG(ss);

    %% ---- PSI ----
    statB_nt = statB; statB_nt.label{1} = 'brain';
    statS_nt = statS; statS_nt.label{1} = 'spine';
    statEMG_nt = statEMG; statEMG_nt.label{1} = 'EMG';
    alldat_nt = ft_appenddata([], statB_nt, statS_nt, statEMG_nt);

    cfgnt = []; cfgnt.output = 'fourier'; cfgnt.method = 'mtmfft';
    cfgnt.foilim = [2 75]; cfgnt.taper = 'hanning'; cfgnt.keeptrials = 'yes';
    freq_nt = ft_freqanalysis(cfgnt, alldat_nt);

    cfgpsi = []; cfgpsi.method = 'psi';
    cfgpsi.bandwidth = fband(2) - fband(1);
    psi_all = ft_connectivityanalysis(cfgpsi, freq_nt);

    iB_psi  = find(strcmp(psi_all.label,'brain'));
    iS_psi  = find(strcmp(psi_all.label,'spine'));
    iE_psi  = find(strcmp(psi_all.label,'EMG'));
    fbmask  = psi_all.freq >= fband(1) & psi_all.freq <= fband(2);

    results(ss).psi_brainEMG   = sum(squeeze(psi_all.psispctrm(iB_psi, iE_psi, fbmask)));
    results(ss).psi_brainSpine = sum(squeeze(psi_all.psispctrm(iB_psi, iS_psi, fbmask)));
    results(ss).psi_spineEMG   = sum(squeeze(psi_all.psispctrm(iS_psi, iE_psi, fbmask)));

    %% ---- Peak cross-correlation latencies ----
    [results(ss).brainEMG.peak_rho,   results(ss).brainEMG.peak_latency]   = ...
        select_peak(t1_nt, lat_min, lat_max);
    [results(ss).brainSpine.peak_rho, results(ss).brainSpine.peak_latency] = ...
        select_peak(t2_nt, lat_min, lat_max);
    [results(ss).spineEMG.peak_rho,   results(ss).spineEMG.peak_latency]   = ...
        select_peak(t3_nt, lat_min, lat_max);

    results(ss).brainEMG.top_lats   = get_top_peaks(t1_nt, lat_min, lat_max, 3);
    results(ss).brainSpine.top_lats = get_top_peaks(t2_nt, lat_min, lat_max, 3);
    results(ss).spineEMG.top_lats   = get_top_peaks(t3_nt, lat_min, lat_max, 3);

end  % subject loop


%% ---- GMM clustering of peak frequencies ---------------------------------
p1_idx = 1;   % OP00212 is first subject
% Pool peak frequencies across all three pairs to find natural clusters.
% Uses MATLAB fitgmdist; BIC selects number of components.

all_peaks_pooled = [peak_freq_brainEMG; peak_freq_brainSpine; peak_freq_spineEMG];
all_peaks_pooled = all_peaks_pooled(~isnan(all_peaks_pooled));

fprintf('\n=== GMM clustering of peak frequencies ===\n');
fprintf('  Pooled peaks (n=%d): ', numel(all_peaks_pooled));
fprintf('%.1f  ', sort(all_peaks_pooled)');
fprintf('\n');

max_k  = 4;
bic    = nan(1, max_k);
gm_fit = cell(1, max_k);
for k = 1:max_k
    try
        opts = statset('MaxIter',500);
        gm_fit{k} = fitgmdist(all_peaks_pooled, k, ...
            'CovarianceType','full','SharedCovariance',false, ...
            'Replicates',20,'Options',opts);
        bic(k) = gm_fit{k}.BIC;
    catch
        bic(k) = inf;
    end
end
[~, best_k] = min(bic);
gm_best     = gm_fit{best_k};
[~, sort_idx]   = sort(gm_best.mu);
cluster_means   = gm_best.mu(sort_idx);

fprintf('  BIC: '); fprintf('k=%d: %.1f  ',[1:max_k; bic]); fprintf('\n');
fprintf('  Best k = %d\n', best_k);
fprintf('  Cluster means (Hz): '); fprintf('%.1f  ', cluster_means); fprintf('\n');

% Boundaries: midpoints between adjacent sorted means + fband edges
boundaries = [fband(1), ...
    arrayfun(@(i) mean(cluster_means(i:i+1)), 1:best_k-1), ...
    fband(2)];
fprintf('  Boundaries (Hz): '); fprintf('%.1f  ', boundaries); fprintf('\n');

nClusters    = best_k;
cluster_cols = lines(nClusters);

% Fallback to gap method if k=1
if nClusters < 2
    warning('GMM k=1; falling back to largest-gap split.');
    sp_ = sort(all_peaks_pooled); [~,gi] = max(diff(sp_));
    boundaries    = [fband(1), mean(sp_(gi:gi+1)), fband(2)];
    nClusters     = 2;
    cluster_means = [mean(sp_(sp_<boundaries(2))), mean(sp_(sp_>=boundaries(2)))];
    cluster_cols  = lines(2);
end
fprintf('======================================\n\n');

%% ---- Sub-band peak frequencies per pair (GMM clusters) -----------------
peak_freq_clust = nan(nSubs, 3, nClusters);   % subs x pairs x clusters
for ss = 1:nSubs
    fax = results(ss).freq_axis_mt;
    if isempty(fax), continue; end
    pairs_dir = {results(ss).dir_BE, results(ss).dir_BS, results(ss).dir_SE};
    for p = 1:3
        d = pairs_dir{p};
        if isempty(d) || all(isnan(d)), continue; end
        for c = 1:nClusters
            mask_c = fax >= boundaries(c) & fax <= boundaries(c+1);
            f_c    = fax(mask_c);
            if sum(mask_c) < 2, continue; end
            [~, ic] = max(d(mask_c));
            peak_freq_clust(ss,p,c) = f_c(ic);
        end
    end
end

% Legacy aliases used by aligned plot downstream
peak_freq_low  = peak_freq_clust(:,:,1);
peak_freq_high = peak_freq_clust(:,:,nClusters);
beta_split_hz  = boundaries(2);
fband_low      = [boundaries(1)          boundaries(2)];
fband_high     = [boundaries(nClusters)  boundaries(nClusters+1)];

pair_names_r = {'Brain-EMG','Brain-Spine','Spine-EMG'};
fprintf('  Per-pair peak frequencies by cluster:\n');
for c = 1:nClusters
    fprintf('  Cluster %d (%.1f-%.1f Hz, mean=%.1f Hz):\n', ...
        c, boundaries(c), boundaries(c+1), cluster_means(c));
    for p = 1:3
        pf = peak_freq_clust(:,p,c); pf = pf(~isnan(pf));
        if isempty(pf)
            fprintf('    %s: no data\n', pair_names_r{p}); continue
        end
        fprintf('    %s: median=%.1f  range=[%.1f %.1f]  n=%d\n', ...
            pair_names_r{p}, median(pf), min(pf), max(pf), numel(pf));
    end
end

%% ---- Partial coherence reduction by GMM cluster (Brain-EMG | spine) ----
% Option 1: peak full vs partial coherence within each cluster per subject.
% Option 2: group aligned full vs partial spectra (aligned to global BE peak).

% --- Option 1: cluster peak values ---
part_full_clust    = nan(nSubs, nClusters);
part_partial_clust = nan(nSubs, nClusters);

for ss = 1:nSubs
    fax = results(ss).freq_axis_mt;
    if isempty(fax), continue; end
    full_v    = results(ss).full_coh_mt_BE;
    partial_v = results(ss).partial_mt_BE_spine;
    for c = 1:nClusters
        mask_c = fax >= boundaries(c) & fax <= boundaries(c+1);
        if sum(mask_c) < 2, continue; end
        part_full_clust(ss,c)    = max(full_v(mask_c));
        part_partial_clust(ss,c) = max(partial_v(mask_c));
    end
end

part_abs_red  = part_full_clust - part_partial_clust;
part_prop_red = part_abs_red ./ part_full_clust;

fprintf('\n=== Partial coherence reduction by cluster (Brain-EMG | spine, MT) ===\n');
for c = 1:nClusters
    fprintf('  Cluster %d (%.1f-%.1f Hz):\n', c, boundaries(c), boundaries(c+1));
    fprintf('    Sub         Full    Partial  AbsRed  PropRed\n');
    for ss = 1:nSubs
        fprintf('    %-12s %6.4f  %6.4f   %6.4f  %6.3f\n', subs{ss}, ...
            part_full_clust(ss,c), part_partial_clust(ss,c), ...
            part_abs_red(ss,c), part_prop_red(ss,c));
    end
    fprintf('    Group: Full median=%.4f  Partial median=%.4f  PropRed median=%.3f\n', ...
        median(part_full_clust(:,c),'omitnan'), ...
        median(part_partial_clust(:,c),'omitnan'), ...
        median(part_prop_red(:,c),'omitnan'));
end
fprintf('======================================================================\n\n');

% Box plot: full vs partial per cluster
col_full    = col_BE;
col_partial = [0.6 0.2 0.8];

hfig_part_clust = figure('Color','w','Position',[100 100 nClusters*380 380]);
for c = 1:nClusters
    subplot(1,nClusters,c); hold on;
    ttl = sprintf('Cluster %d (%.0f-%.0f Hz)', c, boundaries(c), boundaries(c+1));

    for ss = 1:nSubs
        x = [1 2]; y = [part_full_clust(ss,c), part_partial_clust(ss,c)];
        if any(isnan(y)), continue; end
        lc = [0.7 0.7 0.7]; lw = 0.8;
        if ss == p1_idx, lc = [0 0 0]; lw = 2; end
        plot(x, y*1e3, '-o', 'Color',lc, 'LineWidth',lw, 'MarkerSize',5, ...
            'MarkerFaceColor',lc);
    end
    % group medians
    plot([1 2], [median(part_full_clust(:,c),'omitnan') ...
                 median(part_partial_clust(:,c),'omitnan')]*1e3, ...
        '--s', 'Color',col_BE, 'LineWidth',2.5, 'MarkerSize',9, ...
        'MarkerFaceColor','w');

    set(gca,'XTick',[1 2],'XTickLabel',{'Full','Partial'},'FontSize',13);
    ylabel('Peak coherence (x1e-3)','FontSize',13);
    title(ttl,'FontSize',13); xlim([0.6 2.4]); box off;
end
sgtitle('Brain-EMG coherence: full vs partial (spine removed) by cluster','FontSize',13);
if saveFigs
    savefig(hfig_part_clust, fullfile(fig_dir,'partial_coherence_by_cluster.fig'));
end

% --- Option 2: group aligned full vs partial spectra ---
rel_win_p  = 12;
rel_res_p  = 0.25;
rel_axis_p = (-rel_win_p : rel_res_p : rel_win_p);
nRelF_p    = length(rel_axis_p);

% Align to global BE peak per subject
hfig_part_aligned = figure('Color','w','Position',[100 100 500 380]);
hold on;

full_amat    = nan(nSubs, nRelF_p);
partial_amat = nan(nSubs, nRelF_p);

for ss = 1:nSubs
    fax      = results(ss).freq_axis_mt;
    full_v   = results(ss).full_coh_mt_BE;
    part_v   = results(ss).partial_mt_BE_spine;
    if isempty(fax), continue; end

    % align to global BE peak
    mask_full = fax >= fband(1) & fax <= fband(2);
    [~, ig]   = max(full_v(mask_full));
    fc_s      = fax(mask_full); fc_s = fc_s(ig);

    fax_rel  = fax - fc_s;
    in_range = rel_axis_p >= min(fax_rel) & rel_axis_p <= max(fax_rel);
    if sum(in_range) < 2, continue; end

    full_amat(ss,in_range)    = interp1(fax_rel, full_v,  rel_axis_p(in_range), 'linear');
    partial_amat(ss,in_range) = interp1(fax_rel, part_v,  rel_axis_p(in_range), 'linear');
end

% plot mean full and mean partial with SEM shading
mean_full    = mean(full_amat,    1, 'omitnan');
mean_partial = mean(partial_amat, 1, 'omitnan');
sem_full     = std(full_amat,    0, 1, 'omitnan') / sqrt(nSubs);
sem_partial  = std(partial_amat, 0, 1, 'omitnan') / sqrt(nSubs);
ok_cols      = ~all(isnan(full_amat), 1);

fill([rel_axis_p(ok_cols) fliplr(rel_axis_p(ok_cols))], ...
    [(mean_full(ok_cols)+sem_full(ok_cols)) ...
     fliplr(mean_full(ok_cols)-sem_full(ok_cols))]*1e3, ...
    col_full, 'FaceAlpha',0.2, 'EdgeColor','none');
fill([rel_axis_p(ok_cols) fliplr(rel_axis_p(ok_cols))], ...
    [(mean_partial(ok_cols)+sem_partial(ok_cols)) ...
     fliplr(mean_partial(ok_cols)-sem_partial(ok_cols))]*1e3, ...
    col_partial, 'FaceAlpha',0.2, 'EdgeColor','none');

hf = plot(rel_axis_p, mean_full*1e3,    '-',  'Color',col_full,    'LineWidth',2.5);
hp = plot(rel_axis_p, mean_partial*1e3, '--', 'Color',col_partial, 'LineWidth',2.5);
% cluster boundary lines on x axis
for c = 2:nClusters
    xline(boundaries(c) - median(peak_freq_brainEMG,'omitnan'), ...
        '--','Color',[0.5 0.5 0.5],'Alpha',0.4,'HandleVisibility','off');
end

xline(0,'--k','LineWidth',1.2,'HandleVisibility','off');
xlim([-rel_win_p rel_win_p]);
xlabel('Hz from Brain-EMG peak','FontSize',13);
ylabel('Coherence (x1e-3)','FontSize',13);
hf = plot(nan,nan,'-','Color',col_full,'LineWidth',2);
hp = plot(nan,nan,'--','Color',col_partial,'LineWidth',2);
legend([hf hp],{'Full (BE)','Partial (spine removed)'},'Location','northeast','Box','off','FontSize',11);
title('Brain-EMG: full vs partial coherence (aligned)','FontSize',13);
set(gca,'FontSize',13); box off;

if saveFigs
    savefig(hfig_part_aligned, fullfile(fig_dir,'partial_coherence_aligned.fig'));
end
pair_names   = {'Brain-EMG','Brain-Spine','Spine-EMG'};
pair_cols_bp = {col_BE, col_BS, col_SE};

hfig_bp = figure('Color','w','Position',[100 100 nClusters*380 420]);
for c = 1:nClusters
    subplot(1,nClusters,c); hold on;
    pf_mat = peak_freq_clust(:,:,c);
    ttl    = sprintf('Cluster %d (%.0f-%.0f Hz)', c, boundaries(c), boundaries(c+1));

    for p = 1:3
        vals = pf_mat(:,p); ok = ~isnan(vals);
        col_p = pair_cols_bp{p};
        if sum(ok) < 2, continue; end
        q25=prctile(vals(ok),25); q75=prctile(vals(ok),75);
        med=median(vals(ok));
        q05=prctile(vals(ok),5);  q95=prctile(vals(ok),95);
        plot([p p],[q05 q25],'-','Color',col_p,'LineWidth',1.2);
        plot([p p],[q75 q95],'-','Color',col_p,'LineWidth',1.2);
        rectangle('Position',[p-0.2 q25 0.4 q75-q25], ...
            'EdgeColor',col_p,'LineWidth',1.5,'FaceColor',[col_p 0.15]);
        plot([p-0.2 p+0.2],[med med],'-','Color',col_p,'LineWidth',2.5);
        non_p1 = find(ok & (1:nSubs)'~=p1_idx);
        jitter = (rand(sum(ok),1)-0.5)*0.15;
        scatter(p+jitter(ok & (1:nSubs)'~=p1_idx), vals(non_p1), ...
            30,col_p,'filled','MarkerFaceAlpha',0.5);
        if ok(p1_idx)
            scatter(p, vals(p1_idx), 60,'k','filled','Marker','d');
        end
    end
    for ss = 1:nSubs
        row = pf_mat(ss,:); if any(isnan(row)), continue; end
        lc=[0.7 0.7 0.7]; lw=0.8;
        if ss==p1_idx, lc=[0 0 0]; lw=1.5; end
        plot(1:3,row,'-','Color',lc,'LineWidth',lw);
    end
    set(gca,'XTick',1:3,'XTickLabel',pair_names,'FontSize',13);
    ylabel('Peak frequency (Hz)','FontSize',13);
    title(ttl,'FontSize',13); box off;
    ylim([boundaries(c)-1 boundaries(c+1)+1]);
end
sgtitle('Peak frequency tracking across signal pairs (GMM clusters)','FontSize',13);
if saveFigs
    savefig(hfig_bp, fullfile(fig_dir,'group_subband_peak_frequencies.fig'));
end

%% ---- KDE + scatter: GMM cluster version ---------------------------------
kde_bw = 1.5;
f_eval = linspace(fband(1)-2, fband(2)+2, 300);

hfig_clust = figure('Color','w','Position',[100 100 1100 600]);

% Row 1: KDE per pair
for p = 1:3
    subplot(2,4,p); hold on;
    max_kde = 0;
    for c = 1:nClusters
        pf_c = peak_freq_clust(:,p,c); ok_c = ~isnan(pf_c);
        if sum(ok_c) < 2, continue; end
        kde_c = zeros(size(f_eval));
        for ss = 1:nSubs
            if ok_c(ss), kde_c = kde_c + normpdf(f_eval, pf_c(ss), kde_bw); end
        end
        kde_c   = kde_c / sum(ok_c);
        max_kde = max(max_kde, max(kde_c));
        fill([f_eval fliplr(f_eval)],[kde_c zeros(size(kde_c))], ...
            cluster_cols(c,:),'FaceAlpha',0.25,'EdgeColor','none');
        plot(f_eval, kde_c,'Color',cluster_cols(c,:),'LineWidth',2);
    end
    rug_y = max_kde * 0.05;
    for c = 1:nClusters
        pf_c = peak_freq_clust(:,p,c); ok_c = ~isnan(pf_c);
        for ss = 1:nSubs
            if ~ok_c(ss), continue; end
            lw_r = 0.8; if ss==p1_idx, lw_r=2; end
            plot([pf_c(ss) pf_c(ss)],[0 rug_y],'-', ...
                'Color',cluster_cols(c,:),'LineWidth',lw_r);
        end
    end
    for bi = 2:nClusters
        xline(boundaries(bi),'--k','Alpha',0.3,'HandleVisibility','off');
    end
    xlim([fband(1)-1 fband(2)+1]);
    xlabel('Peak frequency (Hz)','FontSize',11);
    ylabel('Density','FontSize',11);
    title(pair_names_r{p},'FontSize',12);
    set(gca,'FontSize',11); box off;
    if p == 1
        leg_e = arrayfun(@(c) sprintf('Cluster %d (%.0f-%.0f Hz)', ...
            c,boundaries(c),boundaries(c+1)), 1:nClusters,'UniformOutput',false);
        legend(leg_e,'Location','northwest','Box','off','FontSize',9);
    end
end

% Row 2: scatter BE vs BS, BE vs SE, BS vs SE
markers = {'o','d','s','^'};
scatter_cfg = {[1 2],'Brain-EMG vs Brain-Spine','Brain-EMG peak (Hz)','Brain-Spine peak (Hz)'; ...
               [1 3],'Brain-EMG vs Spine-EMG',  'Brain-EMG peak (Hz)','Spine-EMG peak (Hz)'; ...
               [2 3],'Brain-Spine vs Spine-EMG','Brain-Spine peak (Hz)','Spine-EMG peak (Hz)'};
for sp = 1:3
    subplot(2,4,4+sp); hold on;
    pi1 = scatter_cfg{sp,1}(1); pi2 = scatter_cfg{sp,1}(2);
    dlim = [fband(1)-1 fband(2)+1];
    plot(dlim,dlim,'--','Color',[0.7 0.7 0.7],'LineWidth',1,'HandleVisibility','off');
    for bi = 2:nClusters
        xline(boundaries(bi),'--k','Alpha',0.2,'HandleVisibility','off');
        yline(boundaries(bi),'--k','Alpha',0.2,'HandleVisibility','off');
    end
    for ss = 1:nSubs
        ms=60; lw=1; if ss==p1_idx, ms=100; lw=2; end
        xv = squeeze(peak_freq_clust(ss,pi1,:));
        yv = squeeze(peak_freq_clust(ss,pi2,:));
        ok_xy = ~isnan(xv) & ~isnan(yv);
        if sum(ok_xy)>1
            plot(xv(ok_xy),yv(ok_xy),'-','Color',[0.8 0.8 0.8], ...
                'LineWidth',0.8,'HandleVisibility','off');
        end
        for c = 1:nClusters
            xc=peak_freq_clust(ss,pi1,c); yc=peak_freq_clust(ss,pi2,c);
            if isnan(xc)||isnan(yc), continue; end
            scatter(xc,yc,ms,cluster_cols(c,:),markers{c},'filled', ...
                'MarkerEdgeColor','w','LineWidth',lw, ...
                'MarkerFaceAlpha',0.75,'HandleVisibility','off');
        end
    end
    hleg = gobjects(nClusters,1);
    for c = 1:nClusters
        hleg(c) = scatter(nan,nan,60,cluster_cols(c,:),markers{c},'filled', ...
            'MarkerFaceAlpha',0.75);
    end
    legend(hleg, arrayfun(@(c) sprintf('Cluster %d',c), ...
        1:nClusters,'UniformOutput',false), ...
        'Location','northwest','Box','off','FontSize',10);
    xlim(dlim); ylim(dlim); axis square;
    xlabel(scatter_cfg{sp,3},'FontSize',11);
    ylabel(scatter_cfg{sp,4},'FontSize',11);
    title(scatter_cfg{sp,2},'FontSize',12);
    set(gca,'FontSize',11); box off;
end

subplot(2,4,8); axis off;
ann = [{'GMM clusters (BIC-selected):'}, ...
    arrayfun(@(c) sprintf('C%d: %.0f-%.0f Hz (mean %.1f)', ...
        c,boundaries(c),boundaries(c+1),cluster_means(c)), ...
        1:nClusters,'UniformOutput',false), ...
    {sprintf('KDE bw=%.1f Hz | k=%d | Large=P1', kde_bw, nClusters)}];
text(0.5,0.55,ann,'HorizontalAlignment','center','FontSize',10,'Units','normalized');

sgtitle('Peak frequency clustering (GMM)','FontSize',13);
if saveFigs
    savefig(hfig_clust, fullfile(fig_dir,'group_peak_frequency_clustering.fig'));
end

%% ---- Group: individual traces aligned to global beta peak ---------------
% One panel per pair. Each subject aligned to their own global peak in
% 10-35 Hz. Individual traces only; P1 in black. No group average line.

rel_win  = 12;
rel_res  = 0.25;
rel_axis = (-rel_win : rel_res : rel_win);
nRelF    = length(rel_axis);

pair_data_g  = {all_dir_brainEMG_stat, all_dir_brainSpine_stat, all_dir_spineEMG_stat};
pair_cols_g  = {col_BE, col_BS, col_SE};
pair_names_g = {'Brain-EMG', 'Brain-Spine', 'Spine-EMG'};

% Global peak per subject per pair (max in full beta band)
global_peak_freq = nan(nSubs, 3);
for ss = 1:nSubs
    fax = results(ss).freq_axis_mt;
    if isempty(fax), continue; end
    mask_full = fax >= fband(1) & fax <= fband(2);
    f_full    = fax(mask_full);
    pairs_dir = {results(ss).dir_BE, results(ss).dir_BS, results(ss).dir_SE};
    for p = 1:3
        d = pairs_dir{p};
        if isempty(d) || all(isnan(d)), continue; end
        [~, ig] = max(d(mask_full));
        global_peak_freq(ss,p) = f_full(ig);
    end
end

hfig_aligned = figure('Color','w','Position',[100 100 900 320]);

for p = 1:3
    subplot(1,3,p); hold on;
    coh_mat     = pair_data_g{p};
    col_p       = pair_cols_g{p};
    align_freqs = global_peak_freq(:,p);

    for ss = 1:nSubs
        fc_s = align_freqs(ss);
        if isnan(fc_s), continue; end
        fax_rel  = freq_axis - fc_s;
        in_range = rel_axis >= min(fax_rel) & rel_axis <= max(fax_rel);
        if sum(in_range) < 2, continue; end
        y = interp1(fax_rel, coh_mat(ss,:), rel_axis(in_range), 'linear') * 1e3;

        if ss == p1_idx
            plot(rel_axis(in_range), y, 'Color','k', 'LineWidth',2);
        else
            plot(rel_axis(in_range), y, 'Color',[col_p 0.5], 'LineWidth',1.2);
        end
    end

    xline(0,'--k','LineWidth',1.2,'HandleVisibility','off');
    xlim([-rel_win rel_win]);
    xlabel('Hz from peak','FontSize',12);
    ylabel('Coherence (x1e-3)','FontSize',12);
    title(pair_names_g{p},'FontSize',13);
    set(gca,'FontSize',12); box off;
end

sgtitle('Group coherence spectra aligned to individual beta peak','FontSize',13);
if saveFigs
    savefig(hfig_aligned, fullfile(fig_dir,'group_coherence_spectra_aligned.fig'));
end
okB = isfinite(Pstat_brain) & isfinite(Prest_brain);
[~,pb,~,sb] = ttest(log(Pstat_brain(okB)), log(Prest_brain(okB)));
fprintf('  Brain: n=%d, t(%d)=%.3f, p=%.4g\n', sum(okB), sb.df, sb.tstat, pb);
okS = isfinite(Pstat_spine) & isfinite(Prest_spine);
[~,ps,~,ss_] = ttest(log(Pstat_spine(okS)), log(Prest_spine(okS)));
fprintf('  Spine: n=%d, t(%d)=%.3f, p=%.4g\n', sum(okS), ss_.df, ss_.tstat, ps);

fprintf('  --- Periodic power (FOOOF peaks, 10-35 Hz) ---\n');
okBf = isfinite(Pstat_brain_fooof) & isfinite(Prest_brain_fooof);
okSf = isfinite(Pstat_spine_fooof) & isfinite(Prest_spine_fooof);
if sum(okBf) > 1
    [~,pbf,~,sbf] = ttest(Pstat_brain_fooof(okBf), Prest_brain_fooof(okBf));
    fprintf('  Brain (periodic): n=%d, t(%d)=%.3f, p=%.4g\n', sum(okBf), sbf.df, sbf.tstat, pbf);
    fprintf('    Contraction: median=%.4f IQR=[%.4f %.4f]\n', ...
        median(Pstat_brain_fooof(okBf)), prctile(Pstat_brain_fooof(okBf),25), prctile(Pstat_brain_fooof(okBf),75));
    fprintf('    Rest:        median=%.4f IQR=[%.4f %.4f]\n', ...
        median(Prest_brain_fooof(okBf)), prctile(Prest_brain_fooof(okBf),25), prctile(Prest_brain_fooof(okBf),75));
else
    fprintf('  Brain (periodic): insufficient data (FOOOF may not have run)\n');
end
if sum(okSf) > 1
    [~,psf,~,ssf] = ttest(Pstat_spine_fooof(okSf), Prest_spine_fooof(okSf));
    fprintf('  Spine (periodic): n=%d, t(%d)=%.3f, p=%.4g\n', sum(okSf), ssf.df, ssf.tstat, psf);
    fprintf('    Contraction: median=%.4f IQR=[%.4f %.4f]\n', ...
        median(Pstat_spine_fooof(okSf)), prctile(Pstat_spine_fooof(okSf),25), prctile(Pstat_spine_fooof(okSf),75));
    fprintf('    Rest:        median=%.4f IQR=[%.4f %.4f]\n', ...
        median(Prest_spine_fooof(okSf)), prctile(Prest_spine_fooof(okSf),25), prctile(Prest_spine_fooof(okSf),75));
else
    fprintf('  Spine (periodic): insufficient data (FOOOF may not have run)\n');
end

%% ---- Export supplementary tables: peak frequency and coherence ----------

% --- Table 1a: Peak frequency (Hz, MT) ---
peak_freq_table_full = array2table(...
    [peak_freq_brainEMG, peak_freq_brainSpine, peak_freq_spineEMG], ...
    'VariableNames', {'BrainEMG_Hz','BrainSpine_Hz','SpineEMG_Hz'}, ...
    'RowNames', subs);

freq_summary = [...
    median(peak_freq_brainEMG,'omitnan'), ...
    median(peak_freq_brainSpine,'omitnan'), ...
    median(peak_freq_spineEMG,'omitnan'); ...
    mad(peak_freq_brainEMG,1), ...
    mad(peak_freq_brainSpine,1), ...
    mad(peak_freq_spineEMG,1)];

freq_summary_table = array2table(freq_summary, ...
    'VariableNames', {'BrainEMG_Hz','BrainSpine_Hz','SpineEMG_Hz'}, ...
    'RowNames', {'Median','MAD'});

supp_table_1a = [peak_freq_table_full; freq_summary_table];

fprintf('\n  Supplementary Table 1a: Peak coherence frequency (Hz)\n');
disp(supp_table_1a);

csv_path_1a = fullfile(fig_dir, 'SuppTable1a_peak_frequency.csv');
try
    writetable(supp_table_1a, csv_path_1a, 'WriteRowNames', true);
    fprintf('  Saved: %s\n', csv_path_1a);
catch ME
    warning('Could not save Table 1a (file may be open in Excel): %s', ME.message);
end

% --- Table 1b: Peak coherence value (R2, MT) ---
peak_coh_table_full = array2table(...
    [peak_coh_brainEMG, peak_coh_brainSpine, peak_coh_spineEMG], ...
    'VariableNames', {'BrainEMG_R2','BrainSpine_R2','SpineEMG_R2'}, ...
    'RowNames', subs);

coh_summary = [...
    median(peak_coh_brainEMG,'omitnan'), ...
    median(peak_coh_brainSpine,'omitnan'), ...
    median(peak_coh_spineEMG,'omitnan'); ...
    mad(peak_coh_brainEMG,1), ...
    mad(peak_coh_brainSpine,1), ...
    mad(peak_coh_spineEMG,1)];

coh_summary_table = array2table(coh_summary, ...
    'VariableNames', {'BrainEMG_R2','BrainSpine_R2','SpineEMG_R2'}, ...
    'RowNames', {'Median','MAD'});

supp_table_1b = [peak_coh_table_full; coh_summary_table];

fprintf('\n  Supplementary Table 1b: Peak coherence value (R²)\n');
disp(supp_table_1b);

csv_path_1b = fullfile(fig_dir, 'SuppTable1b_peak_coherence.csv');
try
    writetable(supp_table_1b, csv_path_1b, 'WriteRowNames', true);
    fprintf('  Saved: %s\n', csv_path_1b);
catch ME
    warning('Could not save Table 1b (file may be open in Excel): %s', ME.message);
end
% --- Statistical tests ---
% Friedman test: non-parametric repeated measures, tests whether peak
% frequency differs across the three signal pairs (within-subject factor).
% Each subject is a block; the three pairs are the treatments.
% Requires complete data (no NaN rows).
freq_mat = [peak_freq_brainEMG, peak_freq_brainSpine, peak_freq_spineEMG];
complete = all(isfinite(freq_mat), 2);
freq_mat_complete = freq_mat(complete, :);
subs_complete     = subs(complete);
n_complete        = sum(complete);

fprintf('\n=== Peak frequency statistics ===\n');
fprintf('  N subjects with complete data across all 3 pairs: %d\n', n_complete);

% Friedman test across pairs
if n_complete >= 3
    [p_fried, tbl_fried, stats_fried] = friedman(freq_mat_complete, 1, 'off');
    fprintf('\n  Friedman test (pairs as within-subject factor):\n');
    fprintf('    chi2(%d) = %.3f, p = %.4g\n', ...
        tbl_fried{2,3}, tbl_fried{2,2}, p_fried);

    % Post-hoc pairwise Wilcoxon signed-rank tests if significant
    pair_names = {'Brain-EMG', 'Brain-Spine', 'Spine-EMG'};
    combos = nchoosek(1:3, 2);
    fprintf('\n  Post-hoc pairwise Wilcoxon signed-rank tests (uncorrected):\n');
    fprintf('  %-20s %-20s   p\n', 'Pair A', 'Pair B');
    p_posthoc = nan(size(combos,1),1);
    for c = 1:size(combos,1)
        a = combos(c,1); b = combos(c,2);
        [p_wil, ~] = signrank(freq_mat_complete(:,a), freq_mat_complete(:,b));
        p_posthoc(c) = p_wil;
        fprintf('  %-20s %-20s   %.4g\n', pair_names{a}, pair_names{b}, p_wil);
    end

    % Bonferroni correction
    fprintf('\n  Post-hoc Bonferroni-corrected alpha (3 comparisons): %.4f\n', 0.05/3);
    fprintf('  Significant after correction:\n');
    for c = 1:size(combos,1)
        a = combos(c,1); b = combos(c,2);
        sig = p_posthoc(c) < (0.05/3);
        fprintf('    %s vs %s: p=%.4g [%s]\n', ...
            pair_names{a}, pair_names{b}, p_posthoc(c), ...
            conditional_str(sig, 'significant', 'ns'));
    end
else
    fprintf('  Insufficient complete cases for Friedman test (need >= 3).\n');
end

% Descriptive summary per pair
fprintf('\n  Descriptive summary (Hz):\n');
fprintf('  %-15s  Median   MAD    Min    Max\n', 'Pair');
for pp = 1:3
    v = freq_mat(:,pp);
    fprintf('  %-15s  %6.2f  %5.2f  %5.2f  %5.2f\n', ...
        pair_names{pp}, median(v,'omitnan'), mad(v,1), ...
        min(v,[],'omitnan'), max(v,[],'omitnan'));
end
fprintf('================================================\n');

%% ---- Collect PSI --------------------------------------------------------
for ss = 1:nSubs
    if isfield(results(ss),'psi_brainEMG')
        psi_brainEMG(ss)   = results(ss).psi_brainEMG;
        psi_brainSpine(ss) = results(ss).psi_brainSpine;
        psi_spineEMG(ss)   = results(ss).psi_spineEMG;
    end
end


%% ---- Group: directed coherence bar plot ---------------------------------
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
ylabel('Coherence area (10-35 Hz)','FontSize',14);
ylim([0 max(data_bar(:))*1.3]);
legend({'Forward','Reverse'},'Location','northwest','FontSize',14); box off;
title('Group directed coherence', 'Interpreter','none','FontSize',14);
if saveFigs
    savefig(hfig_bar, fullfile(fig_dir, 'group_directed_coherence.fig'));
end


%% ---- Group: directionality comparison (Halliday R2 vs PSI) --------------
pair_labels_dir = {'Brain-EMG','Brain-Spine','Spine-EMG'};
pair_cols_dir   = {col_BE, col_BS, col_SE};

hal_fwd = nan(nSubs,3);
hal_rev = nan(nSubs,3);
for ss = 1:nSubs
    hal_fwd(ss,:) = [results(ss).brainEMG.forward_area, ...
                     results(ss).brainSpine.forward_area, ...
                     results(ss).spineEMG.forward_area];
    hal_rev(ss,:) = [results(ss).brainEMG.reverse_area, ...
                     results(ss).brainSpine.reverse_area, ...
                     results(ss).spineEMG.reverse_area];
end

psi_mat = [psi_brainEMG, psi_brainSpine, psi_spineEMG];

hfig_dir_comp = figure('Color','w','Position',[100 100 900 560]);

for pp = 1:3
    col_p = pair_cols_dir{pp};
    jit   = 0.08;

    subplot(2,3,pp); hold on;
    for ss = 1:nSubs
        xj = [1+jit*randn, 2+jit*randn];
        yj = [hal_fwd(ss,pp), hal_rev(ss,pp)];
        if ss == p1
            plot(xj, yj, '-o', 'Color',[0.1 0.1 0.1],'LineWidth',2.5,...
                'MarkerSize',8,'MarkerFaceColor',[0.1 0.1 0.1],'MarkerEdgeColor','w');
        else
            plot(xj, yj, '-','Color',[0.75 0.75 0.75],'LineWidth',0.8);
            scatter(xj(1),yj(1),35,[0.2 0.4 0.8],'filled','MarkerFaceAlpha',0.55,'MarkerEdgeColor','none');
            scatter(xj(2),yj(2),35,[0.9 0.3 0.3],'filled','MarkerFaceAlpha',0.55,'MarkerEdgeColor','none');
        end
    end
    plot(1, median(hal_fwd(:,pp),'omitnan'), 's','MarkerSize',12,...
    'MarkerFaceColor',[0.2 0.4 0.8],'MarkerEdgeColor','k','LineWidth',1.5);
    plot(2, median(hal_rev(:,pp),'omitnan'), 's','MarkerSize',12,...
    'MarkerFaceColor',[0.9 0.3 0.3],'MarkerEdgeColor','k','LineWidth',1.5);
    set(gca,'XTick',[1 2],'XTickLabel',{'Forward','Reverse'},'FontSize',14);
    ylabel('Coherence area','FontSize',14);
    title(sprintf('%s\nHalliday R2', pair_labels_dir{pp}),'FontSize',14);
    ylim([0 max([hal_fwd(:,pp);hal_rev(:,pp)])*1.3]); xlim([0.5 2.5]); box off;

    subplot(2,3,3+pp); hold on;
    yline(0,'-k','LineWidth',1,'HandleVisibility','off');
    for ss = 1:nSubs
        xj = 1 + jit*randn;
        yj = psi_mat(ss,pp);
        if isnan(yj), continue; end
        if ss == p1
            plot(xj, yj, 'o','Color',[0.1 0.1 0.1],'LineWidth',2.5,...
                'MarkerSize',10,'MarkerFaceColor',[0.1 0.1 0.1],'MarkerEdgeColor','w');
        else
        scatter(xj, yj, 45, [0.2 0.7 0.3],'filled','MarkerFaceAlpha',0.55,'MarkerEdgeColor','none');
        end
    end
    med_psi = median(psi_mat(:,pp),'omitnan');
    mad_psi = mad(psi_mat(:,pp),1);
    plot(1, med_psi, 's','MarkerSize',12,...
    'MarkerFaceColor',[0.2 0.7 0.3],'MarkerEdgeColor','k','LineWidth',1.5);
    errorbar(1, med_psi, mad_psi,'k','LineWidth',1.5,'CapSize',8);
    set(gca,'XTick',1,'XTickLabel',{'PSI'},'FontSize',14);
    ylabel('Summed PSI (beta band)','FontSize',14);
    title(sprintf('%s\nPSI', pair_labels_dir{pp}),'FontSize',14);
    psi_lim = max(abs(psi_mat(:,pp)),[],'omitnan')*1.4;
    if isnan(psi_lim) || psi_lim == 0, psi_lim = 1; end
    ylim([-psi_lim psi_lim]); xlim([0.5 1.5]); box off;
end

sgtitle('Directionality: Halliday R2 (top) vs PSI (bottom)', 'FontSize',14,'Interpreter','none');
if saveFigs
    savefig(hfig_dir_comp, fullfile(fig_dir, 'group_directionality_comparison.fig'));
end


%% ---- Group: peak latency plot -------------------------------------------
brainEMG_lat   = abs(arrayfun(@(s) s.brainEMG.peak_latency,   results));
brainSpine_lat = abs(arrayfun(@(s) s.brainSpine.peak_latency, results));
spineEMG_lat   = abs(arrayfun(@(s) s.spineEMG.peak_latency,   results));

lightGrey  = [0.75 0.75 0.75];
darkGrey   = [0.25 0.25 0.25];
medianGrey = [0.10 0.10 0.10];
x_lat      = [1 2 3];

hfig_lat = figure('Color','w'); hold on;
hP1  = [];
hMed = [];

for ss = 1:nSubs
    y = [brainEMG_lat(ss), brainSpine_lat(ss), spineEMG_lat(ss)];
    if any(isnan(y)), continue; end
    if ss == p1
        hP1 = plot(x_lat, y, '-o', ...
            'LineWidth',3.5,'MarkerSize',9, ...
            'Color',darkGrey,'MarkerFaceColor',darkGrey,'MarkerEdgeColor','w', ...
            'DisplayName', subs{ss});
    else
        plot(x_lat, y, '-o','Color',lightGrey,'LineWidth',0.8, ...
            'MarkerSize',6,'MarkerFaceColor',lightGrey,'MarkerEdgeColor','w', ...
            'DisplayName', subs{ss});
        text(3.05, y(3), subs{ss}, 'FontSize', 8, 'Color', lightGrey, ...
            'VerticalAlignment','middle');
    end
end

scatter(ones(nSubs,1),   brainEMG_lat',   60, lightGrey, 'filled','MarkerFaceAlpha',0.5,'MarkerEdgeColor','none');
scatter(2*ones(nSubs,1), brainSpine_lat', 60, lightGrey, 'filled','MarkerFaceAlpha',0.5,'MarkerEdgeColor','none');
scatter(3*ones(nSubs,1), spineEMG_lat',   60, lightGrey, 'filled','MarkerFaceAlpha',0.5,'MarkerEdgeColor','none');

med_vals = [median(brainEMG_lat,'omitnan'), ...
            median(brainSpine_lat,'omitnan'), ...
            median(spineEMG_lat,'omitnan')];
mad_vals = [mad(brainEMG_lat,1), mad(brainSpine_lat,1), mad(spineEMG_lat,1)];

hMed = plot(x_lat, med_vals, '--s', ...
    'LineWidth',2.5,'MarkerSize',9, ...
    'Color',medianGrey,'MarkerFaceColor','w','MarkerEdgeColor',medianGrey);
errorbar(x_lat, med_vals, mad_vals, ...
    'Color',medianGrey,'LineStyle','none','LineWidth',1.5,'CapSize',10);

leg_handles = [];
leg_labels  = {};
if ~isempty(hP1)
    leg_handles = [leg_handles, hP1];
    leg_labels  = [leg_labels, {'Participant 1'}];
end
leg_handles = [leg_handles, hMed];
leg_labels  = [leg_labels,  {'Median +/- MAD'}];
legend(leg_handles, leg_labels, 'Location','northwest','Box','off','FontSize',14);

xlim([0.5 3.8]); xticks(x_lat);
xticklabels({'Brain<->EMG','Brain<->Spine','Spine<->EMG'});
ylabel('Peak latency (ms)','FontSize',14); set(gca,'FontSize',14); grid on; box on;
title('Cross-correlation peak latencies', 'Interpreter','none','FontSize',14);
if saveFigs
    savefig(hfig_lat, fullfile(fig_dir, 'group_peak_latencies.fig'));
end


%% ---- Print latency summary ----------------------------------------------
fprintf('\n=== Peak cross-correlation latencies (ms) ===\n');
fprintf('  Pair            P1      Median   MAD\n');
fprintf('  Brain<->EMG   %5.1f   %5.1f    %5.1f\n', ...
    brainEMG_lat(p1), med_vals(1), mad_vals(1));
fprintf('  Brain<->Spine %5.1f   %5.1f    %5.1f\n', ...
    brainSpine_lat(p1), med_vals(2), mad_vals(2));
fprintf('  Spine<->EMG   %5.1f   %5.1f    %5.1f\n', ...
    spineEMG_lat(p1), med_vals(3), mad_vals(3));
fprintf('=======================================================\n');

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
    valid_idx = t(:,1) >= lat_min & t(:,1) <= lat_max;
    rho_vals  = t(valid_idx, 3);
    lag_vals  = t(valid_idx, 1);

    [~, locs] = findpeaks(abs(rho_vals));
    if isempty(locs), locs = (1:numel(rho_vals))'; end
    locs = locs(abs(lag_vals(locs)) > 2);
    if isempty(locs)
        [~, locs] = max(abs(rho_vals));
    end

    [~, ord] = sort(abs(rho_vals(locs)), 'descend');
    locs_sorted = locs(ord);
    nOut = min(nPeaks, numel(locs_sorted));
    lats = lag_vals(locs_sorted(1:nOut));
end

function [peak_rho, peak_latency] = select_peak(t, lat_min, lat_max, zero_excl_ms)
    valid_idx = t(:,1) >= lat_min & t(:,1) <= lat_max;
    rho_vals = t(valid_idx,3);
    lag_vals = t(valid_idx,1);

    if nargin < 4 || isempty(zero_excl_ms)
        zero_excl_ms = 2;
    end

    [~, locs] = findpeaks(abs(rho_vals));

    if isempty(locs)
        locs = (1:numel(rho_vals))';
    end

    locs = locs(abs(lag_vals(locs)) > zero_excl_ms);

    if isempty(locs)
        [~, idx_max] = max(abs(rho_vals));
        peak_rho = rho_vals(idx_max);
        peak_latency = lag_vals(idx_max);
        fprintf('All peaks within +/-%g ms. Returning global max.\n', zero_excl_ms);
        return
    end

    [~, ord] = sort(abs(rho_vals(locs)), 'descend');
    locs_sorted = locs(ord);

    idx_max = locs_sorted(1);
    peak_rho = rho_vals(idx_max);
    peak_latency = lag_vals(idx_max);

    nPrint = min(5, numel(locs_sorted));
    fprintf('\nTop %d peaks outside +/-%g ms:\n', nPrint, zero_excl_ms);
    fprintf('  Rank   Lag (ms)     rho        |rho|\n');
    for k = 1:nPrint
        ii = locs_sorted(k);
        fprintf('   %d     %7.2f    %+8.4g    %8.4g\n', ...
            k, lag_vals(ii), rho_vals(ii), abs(rho_vals(ii)));
    end
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

function [partial_coh, freq_vec, full_coh] = mt_partial_coherence(x, y, z, samp_rate, seg_len, NW)
% MT_PARTIAL_COHERENCE  Multitaper partial coherence R^2_xy.z

K = round(2*NW) - 1;
[tapers, ~] = dpss(seg_len, NW, K);

n_segs  = floor(length(x) / seg_len);
n_freqs = seg_len/2 + 1;

f_xx = zeros(n_freqs, 1);
f_yy = zeros(n_freqs, 1);
f_zz = zeros(n_freqs, 1);
f_xy = zeros(n_freqs, 1);
f_xz = zeros(n_freqs, 1);
f_yz = zeros(n_freqs, 1);

for seg = 1:n_segs
    idx = (seg-1)*seg_len + (1:seg_len);
    xs = x(idx) - mean(x(idx));
    ys = y(idx) - mean(y(idx));
    zs = z(idx) - mean(z(idx));

    for k = 1:K
        tap = tapers(:,k);
        Xk = fft(xs .* tap);  Xk = Xk(1:n_freqs);
        Yk = fft(ys .* tap);  Yk = Yk(1:n_freqs);
        Zk = fft(zs .* tap);  Zk = Zk(1:n_freqs);

        f_xx = f_xx + abs(Xk).^2;
        f_yy = f_yy + abs(Yk).^2;
        f_zz = f_zz + abs(Zk).^2;
        f_xy = f_xy + Xk .* conj(Yk);
        f_xz = f_xz + Xk .* conj(Zk);
        f_yz = f_yz + Yk .* conj(Zk);
    end
end

denom = n_segs * K;
f_xx = f_xx / denom;  f_yy = f_yy / denom;  f_zz = f_zz / denom;
f_xy = f_xy / denom;  f_xz = f_xz / denom;  f_yz = f_yz / denom;

f_xy_z = f_xy - (f_xz .* conj(f_yz)) ./ f_zz;
f_xx_z = f_xx - (abs(f_xz).^2)        ./ f_zz;
f_yy_z = f_yy - (abs(f_yz).^2)        ./ f_zz;

partial_coh = abs(f_xy_z).^2 ./ (f_xx_z .* f_yy_z);
partial_coh = max(0, min(1, real(partial_coh)));

full_coh = abs(f_xy).^2 ./ (f_xx .* f_yy);
full_coh = max(0, min(1, real(full_coh)));

partial_coh(1) = 0;
full_coh(1)    = 0;

freq_vec = (0:n_freqs-1)' * (samp_rate / seg_len);
end
function s = conditional_str(cond, str_true, str_false)
    if cond, s = str_true; else, s = str_false; end
end