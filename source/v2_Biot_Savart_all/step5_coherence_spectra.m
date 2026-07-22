%% step5_brainspine_coherence.m
% Downstream spectra and coherence analysis — BSLaw forward model.
%
% =========================================================================

clear all
close all
clc

%% =========================================================================
%  USER CONFIG
%% =========================================================================
cfg.fieldtrip_path  = 'C:\Users\mspedden\Documents\fieldtrip';
cfg.spm_path        = 'C:\Users\mspedden\Documents\spm';
cfg.bsc_source_path = 'C:\Users\mspedden\Documents\brainspineconnectivity\source';
cfg.neurospec_path  = 'C:\Users\mspedden\Documents\neurospec211NEW\neurospec211';

cfg.data_root = 'C:\spinecoh_data';
cfg.save_dir  = 'C:\Users\mspedden\Documents\brainspine_savetest';

cfg.subs_brain = {'OP00212','OP00213','OP00215','OP00219', ...
    'OP00225','OP00221','OP00224'};
cfg.subs_spine = {'OP00212','OP00213','OP00215','OP00219', ...
    'OP00220','OP00221','OP00224','OP00225','OP00226'};
cfg.subs_coh   = intersect(cfg.subs_brain, cfg.subs_spine,'stable');

ve_configs = {
    'bslaw_prevalence', 'VE_spine_prevalence_sub%s_forspectra_BS.mat', 'VE';
    };

cfg.brain_ve_suffix  = '_brain_pct10';
cfg.spine_ve_suffix  = '_BS';

cfg.fband         = [10 35];   % analysis band (power, coherence, PSI)

cfg.seg_pwr = 11;
cfg.lat_min = -50;
cfg.lat_max =  50;

saveFigs     = 1;
cfg.saveFigs = saveFigs;

%% =========================================================================
%  SETUP
%% =========================================================================
addpath(cfg.bsc_source_path)
addpath(cfg.spm_path)
spm('defaults','EEG')
addpath(cfg.fieldtrip_path)
ft_defaults
addpath(genpath(cfg.neurospec_path))

fprintf('\n=== SPECTRA CONFIG (BS) ===\n')
fprintf('  Coh subs: %d\n', numel(cfg.subs_coh))
fprintf('  Brain VE suffix: %s\n', cfg.brain_ve_suffix)
fprintf('  Spine VE suffix: %s\n', cfg.spine_ve_suffix)
fprintf('  Analysis band:  %.0f-%.0f Hz\n', cfg.fband(1), cfg.fband(2))
fprintf('===========================\n\n')

fprintf('\n>>> STEP B: Coherence spectra\n')
for vi = 1:size(ve_configs,1)
    cfg.spine_ve_pattern = ve_configs{vi,2};
    cfg.spine_ve_varname = ve_configs{vi,3};
    cfg.fig_suffix       = sprintf('_%s', ve_configs{vi,1});
    cfg.fig_dir = fullfile(cfg.save_dir, 'figures', sprintf('spectra_%s', ve_configs{vi,1}));
    if ~exist(cfg.fig_dir,'dir'), mkdir(cfg.fig_dir); end
    fprintf('\n>>> Running VE config: %s\n', ve_configs{vi,1});
    run_coherence_spectra(cfg)
    fprintf('>>> STEP B complete.\n\n')
end
fprintf('\n=== SPECTRA FINISHED ===\n\n')


%% =========================================================================
%  STEP FUNCTIONS
%% =========================================================================

function run_coherence_spectra(cfg)

subs            = cfg.subs_coh;
save_dir        = cfg.save_dir;
data_root       = cfg.data_root;
fband           = cfg.fband;
seg_pwr         = cfg.seg_pwr;
lat_min         = cfg.lat_min;
lat_max         = cfg.lat_max;
saveFigs        = cfg.saveFigs;
fig_dir         = cfg.fig_dir;

nSubs  = numel(subs);
nFreqs = [];
all_coh_brainEMG_stat   = [];
all_coh_brainSpine_stat = [];
all_coh_spineEMG_stat   = [];

peak_coh_brainEMG    = nan(nSubs,1);
peak_coh_brainSpine  = nan(nSubs,1);
peak_coh_spineEMG    = nan(nSubs,1);
peak_coh_brainEMG_nt   = nan(nSubs,1);
peak_coh_brainSpine_nt = nan(nSubs,1);
peak_coh_spineEMG_nt   = nan(nSubs,1);

% Brain-spine peak coherence, full coherence — feeds the
% peak-vs-threshold boxplot at the end of this function.
peak_coh_brainSpine_full = nan(nSubs,1);

freq_axis = [];

Pstat_brain = nan(nSubs,1); Prest_brain = nan(nSubs,1);
Pstat_spine = nan(nSubs,1); Prest_spine = nan(nSubs,1);
Pstat_brain_fooof = nan(nSubs,1); Prest_brain_fooof = nan(nSubs,1);
Pstat_spine_fooof = nan(nSubs,1); Prest_spine_fooof = nan(nSubs,1);

results = struct();
psi_brainEMG   = nan(nSubs,1);
psi_brainSpine = nan(nSubs,1);
psi_spineEMG   = nan(nSubs,1);

col_BE = [0.2 0.4 0.8];
col_BS = [0.8 0.2 0.4];
col_SE = [0.2 0.7 0.3];
col_full    = col_BE;
col_partial = [0.6 0.2 0.8];

M_segments = nan;

% P1 trial-level power storage
p1_brain_stat_trials = [];
p1_brain_rest_trials = [];
p1_spine_stat_trials = [];
p1_spine_rest_trials = [];

%% =========================================================================
%  SUBJECT LOOP
%% =========================================================================
for ss = 1:nSubs
    sub = subs{ss};
    fprintf('  [Step B] Subject %s (%d/%d)\n', sub, ss, nSubs);

    datwithEMGmerged = get_datafile(data_root, sub);
    D     = spm_eeg_load(datwithEMGmerged);
    ftdat = spm2fieldtrip(D);

    %% Load brain VE
    bVE_file = fullfile(save_dir, ...
        sprintf('sub%s_VE_brain_M1%s.mat', sub, cfg.brain_ve_suffix));
    bVE_data = load(bVE_file, 'VE_brain');
    bVE      = bVE_data.VE_brain;

    %% Load spine VE
    sVE_file = fullfile(save_dir, sprintf(cfg.spine_ve_pattern, sub));
    sVE_data = load(sVE_file);
    sVE      = sVE_data.(cfg.spine_ve_varname);

    %% EMG (rectified)
    cfg_ft = []; cfg_ft.channel = 'EXG1';
    EMG = ft_selectdata(cfg_ft, ftdat);
    cfg_ft = []; cfg_ft.rectify = 'yes';
    EMG = ft_preprocessing(cfg_ft, EMG);

    %% Separate conditions
    statidx = find(ftdat.trialinfo==1);
    restidx = find(ftdat.trialinfo==2);
    nTrials = min(length(statidx), length(restidx));

    cfg_ft = []; cfg_ft.trials = statidx(1:nTrials);
    statB   = ft_selectdata(cfg_ft, bVE);
    statS   = ft_selectdata(cfg_ft, sVE);
    statEMG = ft_selectdata(cfg_ft, EMG);

    cfg_ft.trials = restidx(1:nTrials);
    restB   = ft_selectdata(cfg_ft, bVE);
    restS   = ft_selectdata(cfg_ft, sVE);
    restEMG = ft_selectdata(cfg_ft, EMG);

    %% Spectral power
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

    cfgsel = []; cfgsel.frequency = fband; cfgsel.avgoverfreq = 'yes';
    bp_stat = ft_selectdata(cfgsel, freq_stat);
    bp_rest = ft_selectdata(cfgsel, freq_rest);
    m_stat  = squeeze(mean(bp_stat.powspctrm, 1));
    m_rest  = squeeze(mean(bp_rest.powspctrm, 1));
    iB = find(strcmp(bp_stat.label,'brain'));
    iS = find(strcmp(bp_stat.label,'spine'));
    Pstat_brain(ss) = m_stat(iB); Prest_brain(ss) = m_rest(iB);
    Pstat_spine(ss) = m_stat(iS); Prest_spine(ss) = m_rest(iS);

    % P1 trial-level power
    if ss == 1
        p1_brain_stat_trials = squeeze(bp_stat.powspctrm(:, iB));
        p1_brain_rest_trials = squeeze(bp_rest.powspctrm(:, iB));
        p1_spine_stat_trials = squeeze(bp_stat.powspctrm(:, iS));
        p1_spine_rest_trials = squeeze(bp_rest.powspctrm(:, iS));
    end

    %% FOOOF
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
        Pstat_brain_fooof(ss) = m_stat_f(iBf); Prest_brain_fooof(ss) = m_rest_f(iBf);
        Pstat_spine_fooof(ss) = m_stat_f(iSf); Prest_spine_fooof(ss) = m_rest_f(iSf);
    catch ME
        warning('FOOOF failed for %s: %s', sub, ME.message);
    end

    %% NeuroSpec
    statBcont   = [statB.trial{:}];
    statScont   = [statS.trial{:}];
    statEMGcont = abs([statEMG.trial{:}]);
    samp_rate   = ftdat.hdr.Fs;
    opt_str     = 'M4';

    [f1_nt, t1_nt, ~]  = sp2a2_R2_mt(statBcont', statEMGcont', samp_rate, seg_pwr);
    [f2_nt, t2_nt, ~]  = sp2a2_R2_mt(statBcont', statScont',   samp_rate, seg_pwr);
    [f3_nt, t3_nt, ~]  = sp2a2_R2_mt(statScont', statEMGcont', samp_rate, seg_pwr);
    [f1_mt, ~, cl1_mt] = sp2a2_R2_mt(statBcont', statEMGcont', samp_rate, seg_pwr, opt_str);
    [f2_mt, ~, cl2_mt] = sp2a2_R2_mt(statBcont', statScont',   samp_rate, seg_pwr, opt_str);
    [f3_mt, ~, cl3_mt] = sp2a2_R2_mt(statScont', statEMGcont', samp_rate, seg_pwr, opt_str);

    %% P1 coherence spectra — standalone figure
    if strcmp(sub, 'OP00212')
        hfig_p1 = figure('Color','w','Position',[100 100 900 320]);
        fax_p1  = f1_mt(:,1);
        coh_BE  = f1_mt(:,4)*1e3;
        coh_BS  = f2_mt(:,4)*1e3;
        coh_SE  = f3_mt(:,4)*1e3;
        pairs      = {coh_BE, coh_BS, coh_SE};
        cols       = {col_BE, col_BS, col_SE};
        titles_p1  = {'Brain-EMG','Brain-Spine','Spine-EMG'};
        all_coh = [coh_BE; coh_BS; coh_SE];
        ylims   = [0, max(all_coh)*1.15];
        for pp = 1:3
            ax = subplot(1,3,pp); hold on;
            fill([fband(1) fband(2) fband(2) fband(1)], ...
                 [ylims(1) ylims(1) ylims(2) ylims(2)], ...
                 [0.9 0.9 0.9],'FaceAlpha',0.3,'EdgeColor','none');
            plot(fax_p1, pairs{pp}, 'Color',cols{pp}, 'LineWidth',2);
            xlim([2 45]); ylim([0 5]);
            xlabel('Frequency (Hz)','FontSize',12);
            if pp==1, ylabel('Coherence (x1e-3)','FontSize',12); end
            title(titles_p1{pp},'FontSize',12,'Interpreter','none');
            set(ax,'FontSize',11); box off;
        end
        sgtitle('Participant 1 — coherence spectra (BSLaw)','FontSize',13,'Interpreter','none');
        if saveFigs
            savefig(hfig_p1, fullfile(fig_dir,['subOP00212_coherence_spectra' cfg.fig_suffix '.fig']));
            saveas(hfig_p1,  fullfile(fig_dir,['subOP00212_coherence_spectra' cfg.fig_suffix '.png']));
        end
    end

    %% Partial coherence
    mt_NW      = 4;
    mt_seg_len = 2^seg_pwr;
    [fp_mt_BE_spine, ~, fc_mt_BE] = mt_partial_coherence(statBcont', statEMGcont', statScont', samp_rate, mt_seg_len, mt_NW);
    [~, ~, fc_mt_SE] = mt_partial_coherence(statScont', statEMGcont', statBcont',   samp_rate, mt_seg_len, mt_NW);
    [~, ~, fc_mt_BS] = mt_partial_coherence(statBcont', statScont',   statEMGcont', samp_rate, mt_seg_len, mt_NW);

    ns_len = size(f1_mt, 1);
    fp_mt_BE_spine = fp_mt_BE_spine(1:ns_len);
    fc_mt_BE       = fc_mt_BE(1:ns_len);
    fc_mt_SE       = fc_mt_SE(1:ns_len);
    fc_mt_BS       = fc_mt_BS(1:ns_len);

    if isnan(M_segments)
        M_segments = floor(length(statBcont) / 2^seg_pwr);
        fprintf('  NeuroSpec: %d segments\n', M_segments);
    end

    if isempty(nFreqs)
        nFreqs = size(f1_mt, 1);
        all_coh_brainEMG_stat   = nan(nSubs, nFreqs);
        all_coh_brainSpine_stat = nan(nSubs, nFreqs);
        all_coh_spineEMG_stat   = nan(nSubs, nFreqs);
    end

    freq_axis = f1_mt(:,1)';
    coh_BE_full = f1_mt(:,4);
    coh_BS_full = f2_mt(:,4);
    coh_SE_full = f3_mt(:,4);
    all_coh_brainEMG_stat(ss,:)   = coh_BE_full';
    all_coh_brainSpine_stat(ss,:) = coh_BS_full';
    all_coh_spineEMG_stat(ss,:)   = coh_SE_full';

    results(ss).coh_BE = coh_BE_full';
    results(ss).coh_BS = coh_BS_full';
    results(ss).coh_SE = coh_SE_full';

    %% Peak coherence — analysis band (magnitude only, for supplementary table)
    band_mask = freq_axis >= fband(1) & freq_axis <= fband(2);
    peak_coh_brainEMG(ss)   = max(coh_BE_full(band_mask));
    peak_coh_brainSpine(ss) = max(coh_BS_full(band_mask));
    peak_coh_spineEMG(ss)   = max(coh_SE_full(band_mask));

    % Brain-spine peak, full coherence (col 4) — feeds the
    % peak-vs-threshold boxplot below.
    peak_coh_brainSpine_full(ss) = max(f2_mt(band_mask,4));

    %% NT peak coherence (analysis band)
    fax_nt_pk = f1_nt(:,1)';
    band_nt   = fax_nt_pk >= fband(1) & fax_nt_pk <= fband(2);
    coh_BE_nt = f1_nt(:,4);
    coh_BS_nt = f2_nt(:,4);
    coh_SE_nt = f3_nt(:,4);
    peak_coh_brainEMG_nt(ss)   = max(coh_BE_nt(band_nt));
    peak_coh_brainSpine_nt(ss) = max(coh_BS_nt(band_nt));
    peak_coh_spineEMG_nt(ss)   = max(coh_SE_nt(band_nt));

    results(ss).sub                 = sub;
    results(ss).partial_mt_BE_spine = fp_mt_BE_spine';
    results(ss).full_coh_mt_BE      = fc_mt_BE';
    results(ss).full_coh_mt_SE      = fc_mt_SE';
    results(ss).full_coh_mt_BS      = fc_mt_BS';
    results(ss).freq_axis_mt        = f1_mt(:,1)';
    results(ss).conf95_BE           = cl1_mt.ch_c95;
    results(ss).conf95_BS           = cl2_mt.ch_c95;
    results(ss).conf95_SE           = cl3_mt.ch_c95;

    %% Directed coherence areas
    % Always uses forward (col 11) / reverse (col 12) directly and
    % separately — this is a directionality analysis (which way is
    % stronger), not "the" coherence metric, so it is intentionally
    % NOT "the" coherence metric reported elsewhere.
    compute_area = @(fmat, fwd_col, rev_col, band) deal( ...
        trapz(fmat(fmat(:,1)>=band(1) & fmat(:,1)<=band(2),1), ...
        fmat(fmat(:,1)>=band(1) & fmat(:,1)<=band(2),fwd_col)), ...
        trapz(fmat(fmat(:,1)>=band(1) & fmat(:,1)<=band(2),1), ...
        fmat(fmat(:,1)>=band(1) & fmat(:,1)<=band(2),rev_col)) );
    [brainEMG_fwd,  brainEMG_rev]   = compute_area(f1_nt, 11, 12, fband);
    [brainSpine_fwd,brainSpine_rev] = compute_area(f2_nt, 11, 12, fband);
    [spineEMG_fwd,  spineEMG_rev]   = compute_area(f3_nt, 11, 12, fband);

    results(ss).brainEMG.forward_area   = brainEMG_fwd;
    results(ss).brainEMG.reverse_area   = brainEMG_rev;
    results(ss).brainEMG.ratio          = brainEMG_fwd/(brainEMG_fwd+brainEMG_rev);
    results(ss).brainSpine.forward_area = brainSpine_fwd;
    results(ss).brainSpine.reverse_area = brainSpine_rev;
    results(ss).brainSpine.ratio        = brainSpine_fwd/(brainSpine_fwd+brainSpine_rev);
    results(ss).spineEMG.forward_area   = spineEMG_fwd;
    results(ss).spineEMG.reverse_area   = spineEMG_rev;
    results(ss).spineEMG.ratio          = spineEMG_fwd/(spineEMG_fwd+spineEMG_rev);

    %% PSI
    statB_nt   = statB;   statB_nt.label{1}   = 'brain';
    statS_nt   = statS;   statS_nt.label{1}   = 'spine';
    statEMG_nt = statEMG; statEMG_nt.label{1} = 'EMG';
    alldat_nt = ft_appenddata([], statB_nt, statS_nt, statEMG_nt);
    cfgnt = []; cfgnt.output = 'fourier'; cfgnt.method = 'mtmfft';
    cfgnt.foilim = [2 75]; cfgnt.taper = 'hanning'; cfgnt.keeptrials = 'yes';
    freq_nt = ft_freqanalysis(cfgnt, alldat_nt);
    cfgpsi = []; cfgpsi.method = 'psi';
    cfgpsi.bandwidth = fband(2) - fband(1);
    psi_all = ft_connectivityanalysis(cfgpsi, freq_nt);
    iB_psi = find(strcmp(psi_all.label,'brain'));
    iS_psi = find(strcmp(psi_all.label,'spine'));
    iE_psi = find(strcmp(psi_all.label,'EMG'));
    fbmask = psi_all.freq >= fband(1) & psi_all.freq <= fband(2);
    results(ss).psi_brainEMG   = sum(squeeze(psi_all.psispctrm(iB_psi,iE_psi,fbmask)));
    results(ss).psi_brainSpine = sum(squeeze(psi_all.psispctrm(iB_psi,iS_psi,fbmask)));
    results(ss).psi_spineEMG   = sum(squeeze(psi_all.psispctrm(iS_psi,iE_psi,fbmask)));

    %% Peak cross-correlation latencies
    [results(ss).brainEMG.peak_rho,   results(ss).brainEMG.peak_latency]   = select_peak(t1_nt, lat_min, lat_max);
    [results(ss).brainSpine.peak_rho, results(ss).brainSpine.peak_latency] = select_peak(t2_nt, lat_min, lat_max);
    [results(ss).spineEMG.peak_rho,   results(ss).spineEMG.peak_latency]   = select_peak(t3_nt, lat_min, lat_max);
    results(ss).brainEMG.top_lats   = get_top_peaks(t1_nt, lat_min, lat_max, 3);
    results(ss).brainSpine.top_lats = get_top_peaks(t2_nt, lat_min, lat_max, 3);
    results(ss).spineEMG.top_lats   = get_top_peaks(t3_nt, lat_min, lat_max, 3);

end  % subject loop

%% =========================================================================
%  SUPPLEMENTARY: Cord-EMG coherence for spine-only subjects
%  The main loop runs only the 7 brain-cap-spine subjects. Cord-EMG needs no
%  brain VE, so it can also be computed for the spine-only subjects. This
%  brings Cord-EMG to the full n=9 for the raw-spectra grid and the
%  peak-frequency summary; the two brain pairs remain n=7.
%% =========================================================================
extra_subs     = setdiff(cfg.subs_spine, subs, 'stable');
spineEMG_extra = struct('sub', {}, 'coh', {}, 'fax', {}, 'peakf', {}, ...
                        'peak_mt', {}, 'peak_nt', {});
for ee = 1:numel(extra_subs)
    esub = extra_subs{ee};
    fprintf('  [Step B extra] Spine-only subject %s (Cord-EMG only)\n', esub);
    try
        datwithEMGmerged = get_datafile(data_root, esub);
        D     = spm_eeg_load(datwithEMGmerged);
        ftdat = spm2fieldtrip(D);

        sVE_file = fullfile(save_dir, sprintf(cfg.spine_ve_pattern, esub));
        sVE_data = load(sVE_file);
        sVE      = sVE_data.(cfg.spine_ve_varname);

        cfg_ft = []; cfg_ft.channel = 'EXG1';
        EMG = ft_selectdata(cfg_ft, ftdat);
        cfg_ft = []; cfg_ft.rectify = 'yes';
        EMG = ft_preprocessing(cfg_ft, EMG);

        statidx = find(ftdat.trialinfo==1);
        restidx = find(ftdat.trialinfo==2);
        nTrials = min(length(statidx), length(restidx));
        cfg_ft = []; cfg_ft.trials = statidx(1:nTrials);
        statS   = ft_selectdata(cfg_ft, sVE);
        statEMG = ft_selectdata(cfg_ft, EMG);

        statScont   = [statS.trial{:}];
        statEMGcont = abs([statEMG.trial{:}]);
        samp_rate   = ftdat.hdr.Fs;

        [f3_mt, ~, ~] = sp2a2_R2_mt(statScont', statEMGcont', samp_rate, seg_pwr, 'M4');
        fax_e       = f3_mt(:,1)';
        coh_SE_full = f3_mt(:,4)';

        % peak frequency + tapered peak coherence in band (matches main loop)
        band_mask   = fax_e >= fband(1) & fax_e <= fband(2);
        f_band      = fax_e(band_mask);
        [pk_mt, ig] = max(coh_SE_full(band_mask));

        % non-tapered R2 for the NT peak-coherence column (matches main-loop f3_nt)
        [f3_nt, ~, ~] = sp2a2_R2_mt(statScont', statEMGcont', samp_rate, seg_pwr);
        fax_nt  = f3_nt(:,1)';
        band_nt = fax_nt >= fband(1) & fax_nt <= fband(2);
        pk_nt   = max(f3_nt(band_nt,4));

        spineEMG_extra(end+1) = struct('sub', esub, 'coh', coh_SE_full, ...
            'fax', fax_e, 'peakf', f_band(ig), 'peak_mt', pk_mt, 'peak_nt', pk_nt); %#ok<AGROW>
    catch ME
        warning('Cord-EMG failed for spine-only subject %s: %s', esub, ME.message);
        spineEMG_extra(end+1) = struct('sub', esub, 'coh', [], ...
            'fax', [], 'peakf', NaN, 'peak_mt', NaN, 'peak_nt', NaN); %#ok<AGROW>
    end
end

p1_idx = 1;  % Participant 1, highlighted throughout the group figures below

%% =========================================================================
%  FIGURES
%% =========================================================================

%% Partial vs full coherence — paired test across participants
% Same comparison the percent-reduction figure below shows visually (full
% vs partial Brain-EMG coherence, controlling for spine), but reduced to
% one paired value per subject: each subject's own full-coherence peak
% (within fband) vs the partial coherence at that same frequency bin.
% Both come from mt_partial_coherence (same estimator/tapering), so the
% comparison isn't confounded by using two different coherence methods.
% This is a separate technique from sp2a2_R2_mt.
full_peak_BE    = nan(nSubs,1);
partial_at_peak = nan(nSubs,1);
for ss = 1:nSubs
    fax = results(ss).freq_axis_mt;
    if isempty(fax), continue; end
    full_v    = results(ss).full_coh_mt_BE;
    part_v    = results(ss).partial_mt_BE_spine;
    band_mask = fax >= fband(1) & fax <= fband(2);
    f_band    = fax(band_mask);
    full_band = full_v(band_mask);
    part_band = part_v(band_mask);
    [pk_full, ig] = max(full_band);
    full_peak_BE(ss)    = pk_full;
    partial_at_peak(ss) = part_band(ig);
end

ok_part = isfinite(full_peak_BE) & isfinite(partial_at_peak);
n_part  = sum(ok_part);
reduction = full_peak_BE(ok_part) - partial_at_peak(ok_part);

fprintf('\n=== Partial vs full Brain-EMG coherence (controlling for spine) ===\n');
fprintf('  Sub           Full      Partial   Reduction\n');
for ss = 1:nSubs
    if ~ok_part(ss), continue; end
    fprintf('  %-12s  %.4f    %.4f    %.4f\n', subs{ss}, ...
        full_peak_BE(ss), partial_at_peak(ss), full_peak_BE(ss)-partial_at_peak(ss));
end
fprintf('  n=%d, %d/%d subjects show full > partial\n', n_part, sum(reduction>0), n_part);
fprintf('  Median reduction: %.4f (MAD %.4f)\n', median(reduction), mad(reduction,1));

[p_part_signrank, ~, stats_signrank] = signrank(full_peak_BE(ok_part), partial_at_peak(ok_part));
fprintf('  Wilcoxon signed-rank: signedrank=%.1f, p=%.4g\n', stats_signrank.signedrank, p_part_signrank);

[~, p_part_ttest, ~, stats_ttest] = ttest(full_peak_BE(ok_part), partial_at_peak(ok_part));
fprintf('  Paired t-test: t(%d)=%.3f, p=%.4g\n', stats_ttest.df, stats_ttest.tstat, p_part_ttest);
fprintf('=====================================================================\n');

%% Partial vs full coherence — percent reduction (supplementary)
pct_reduction = 100 * (full_peak_BE(ok_part) - partial_at_peak(ok_part)) ./ full_peak_BE(ok_part);
subs_part     = subs(ok_part);

[p_pct_signrank, ~, stats_pct] = signrank(pct_reduction);
med_pct = median(pct_reduction);
mad_pct = mad(pct_reduction, 1);

fprintf('\n=== Partial vs full Brain-EMG coherence — percent reduction ===\n');
fprintf('  Sub           %% reduction\n');
for k = 1:numel(pct_reduction)
    fprintf('  %-12s  %6.1f%%\n', subs_part{k}, pct_reduction(k));
end
fprintf('  Median: %.1f%% (MAD %.1f%%)\n', med_pct, mad_pct);
fprintf('  Wilcoxon signed-rank vs 0: signedrank=%.1f, p=%.4g\n', stats_pct.signedrank, p_pct_signrank);
fprintf('================================================================\n');

try
    pct_table = array2table(pct_reduction(:), 'VariableNames', {'PctReduction'}, 'RowNames', subs_part);
    pct_table = [pct_table; array2table([med_pct; mad_pct], 'VariableNames', {'PctReduction'}, 'RowNames', {'Median','MAD'})];
    writetable(pct_table, fullfile(fig_dir,'SuppTable_partial_coherence_pct_reduction.csv'), 'WriteRowNames', true);
catch ME, warning('Could not save percent-reduction table: %s', ME.message); end

% same per-participant colors as the peak-frequency figure (assumes this
% script's subs list is in the same order as peak_freq_and_shape_BS.m's)
pt_cmap = lines(numel(subs));

jit_pct = 0.08;
hfig_pct = figure('Color','w','Position',[100 100 360 420]); hold on;
yline(0,'-k','LineWidth',1,'HandleVisibility','off');
for k = 1:numel(pct_reduction)
    xj  = 1 + jit_pct*randn;
    idx = find(strcmp(subs, subs_part{k}));
    if strcmp(subs_part{k}, subs{p1_idx})
        plot(xj, pct_reduction(k), 'o','Color','k','LineWidth',2.5,'MarkerSize',12, ...
            'MarkerFaceColor',pt_cmap(idx,:),'MarkerEdgeColor','k');
    else
        scatter(xj, pct_reduction(k), 55, pt_cmap(idx,:), 'filled','MarkerFaceAlpha',0.6,'MarkerEdgeColor','none');
    end
    text(xj + 0.06, pct_reduction(k), sprintf('P%d', idx), ...
        'FontSize', 9, 'Color', pt_cmap(idx,:) * 0.75);
end
plot(1, med_pct, 's','MarkerSize',13,'MarkerFaceColor',col_partial,'MarkerEdgeColor','k','LineWidth',1.5);
errorbar(1, med_pct, mad_pct, 'k','LineWidth',1.5,'CapSize',8);
set(gca,'XTick',1,'XTickLabel',{'Brain-EMG'},'FontSize',13);
ylabel('Coherence reduction after partializing cord (%)','FontSize',12);
title(sprintf('Median %.0f%% reduction (Wilcoxon p=%s)', med_pct, pfmt_short(p_pct_signrank)),'FontSize',12);
xlim([0.5 1.5]); set(gca,'FontSize',12); box off;

if saveFigs
    savefig(hfig_pct, fullfile(fig_dir,['partial_coherence_pct_reduction' cfg.fig_suffix '.fig']));
    saveas(hfig_pct,  fullfile(fig_dir,['partial_coherence_pct_reduction' cfg.fig_suffix '.png']));
end

%% Group coherence spectra aligned to peak
pair_data_g  = {all_coh_brainEMG_stat, all_coh_brainSpine_stat, all_coh_spineEMG_stat};
pair_cols_g  = {col_BE, col_BS, col_SE};
pair_names_g = {'Brain-EMG','Brain-Spine','Spine-EMG'};
rel_win  = 12; rel_res = 0.25;
rel_axis = (-rel_win : rel_res : rel_win);
nRelF    = length(rel_axis);

global_peak_freq = nan(nSubs, 3);
for ss = 1:nSubs
    fax = results(ss).freq_axis_mt;
    if isempty(fax), continue; end
    mask_full = fax >= fband(1) & fax <= fband(2);
    f_full    = fax(mask_full);
    pairs_dir = {results(ss).coh_BE, results(ss).coh_BS, results(ss).coh_SE};
    for p = 1:3
        d = pairs_dir{p};
        if isempty(d) || all(isnan(d)), continue; end
        [~,ig] = max(d(mask_full));
        global_peak_freq(ss,p) = f_full(ig);
    end
end

hfig_aligned = figure('Color','w','Position',[100 100 900 320]);
for p = 1:3
    subplot(1,3,p); hold on;
    coh_mat  = pair_data_g{p};
    col_p    = pair_cols_g{p};
    traces_aligned = nan(nSubs, nRelF);
    for ss = 1:nSubs
        fc_s = global_peak_freq(ss,p);
        if isnan(fc_s), continue; end
        fax_rel  = freq_axis - fc_s;
        in_range = rel_axis >= min(fax_rel) & rel_axis <= max(fax_rel);
        if sum(in_range) < 2, continue; end
        y = interp1(fax_rel, coh_mat(ss,:), rel_axis(in_range), 'linear') * 1e3;
        traces_aligned(ss, in_range) = y;
        if ss == p1_idx, continue; end
        plot(rel_axis(in_range), y, '-','Color',[col_p 0.25],'LineWidth',0.8,'HandleVisibility','off');
    end
    med_trace = median(traces_aligned, 1,'omitnan');
    mad_trace = mad(traces_aligned, 1, 1);
    ok_c      = ~all(isnan(traces_aligned), 1);
    fill([rel_axis(ok_c) fliplr(rel_axis(ok_c))], ...
        [med_trace(ok_c)+mad_trace(ok_c) fliplr(med_trace(ok_c)-mad_trace(ok_c))], ...
        col_p,'FaceAlpha',0.15,'EdgeColor','none','HandleVisibility','off');
    plot(rel_axis(ok_c), med_trace(ok_c), '-','Color',col_p,'LineWidth',2.5,'DisplayName','Median');
    fc_s = global_peak_freq(p1_idx,p);
    if ~isnan(fc_s)
        fax_rel  = freq_axis - fc_s;
        in_range = rel_axis >= min(fax_rel) & rel_axis <= max(fax_rel);
        if sum(in_range) >= 2
            y_p1 = interp1(fax_rel, coh_mat(p1_idx,:), rel_axis(in_range), 'linear') * 1e3;
            plot(rel_axis(in_range), y_p1, 'k-','LineWidth',2,'DisplayName','Participant 1');
        end
    end
    xline(0,'--k','LineWidth',1.2,'HandleVisibility','off');
    xlim([-rel_win rel_win]);
    xlabel('Hz from peak','FontSize',12); ylabel('Coherence (x1e-3)','FontSize',12);
    title(pair_names_g{p},'FontSize',13);
    legend('Location','northeast','Box','off','FontSize',10);
    set(gca,'FontSize',12); box off;
end
sgtitle('Group coherence spectra — aligned to each pair''s own peak','FontSize',13);
if saveFigs
    savefig(hfig_aligned, fullfile(fig_dir,['group_coherence_spectra_global' cfg.fig_suffix '.fig']));
    saveas(hfig_aligned,  fullfile(fig_dir,['group_coherence_spectra_global' cfg.fig_suffix '.png']));
end

%% Raw (unaligned) coherence spectra — N x 3 grid, one row per participant
% Shaped for A4 portrait: rows = participants, columns = pairs (Brain-EMG,
% Brain-Cord, Cord-EMG). The brain columns are empty for the spine-only
% subjects. Coherence line only, restricted to the 10-35 Hz band of interest,
% with each panel scaled individually to its own peak.
%
% Canonical participant order used for BOTH this grid and the peak-frequency
% boxplot below: the 7 brain subjects first, in the SAME order as
% plot_freq_correlation_plots.m and the partial-coherence figure (so P1-P7 =
% the same subjects with the same lines() colours across every figure),
% followed by the 2 spine-only subjects as P8-P9.
grid_extras = setdiff(cfg.subs_spine, subs, 'stable');
grid_subs   = [subs, grid_extras];
nRows       = numel(grid_subs);
xr         = [10 35];  % displayed frequency range = analysis band of interest
pair_names = {'Brain-EMG','Brain-Cord','Cord-EMG'};
pair_cols  = {col_BE, col_BS, col_SE};

grid_coh = cell(nRows, 3);   % {participant, pair} coherence (x1e-3)
grid_fax = cell(nRows, 3);   % matching frequency axis
for rr = 1:nRows
    sub_c = grid_subs{rr};
    bi = find(strcmp(subs, sub_c));
    if ~isempty(bi)
        grid_coh{rr,1} = all_coh_brainEMG_stat(bi,:)   * 1e3; grid_fax{rr,1} = freq_axis;
        grid_coh{rr,2} = all_coh_brainSpine_stat(bi,:) * 1e3; grid_fax{rr,2} = freq_axis;
        grid_coh{rr,3} = all_coh_spineEMG_stat(bi,:)   * 1e3; grid_fax{rr,3} = freq_axis;
    else
        ei = find(strcmp({spineEMG_extra.sub}, sub_c));
        if ~isempty(ei) && ~isempty(spineEMG_extra(ei).coh)
            grid_coh{rr,3} = spineEMG_extra(ei).coh * 1e3; grid_fax{rr,3} = spineEMG_extra(ei).fax;
        end
    end
end

% last populated row per column — carries that column's x tick labels, since
% the bottom row (a spine-only subject) is empty in the brain columns
last_row = zeros(3,1);
for cc = 1:3
    for rr = 1:nRows
        if ~isempty(grid_coh{rr,cc}), last_row(cc) = rr; end
    end
end

hfig_grid = figure('Color','w','Position',[100 30 680 960]);   % A4 portrait
for rr = 1:nRows
    for cc = 1:3
        ax = subplot(nRows, 3, (rr-1)*3 + cc); hold(ax,'on');
        y = grid_coh{rr,cc}; f = grid_fax{rr,cc};
        if isempty(y)
            axis(ax,'off');
        else
            plot(ax, f, y, 'Color', pair_cols{cc}, 'LineWidth', 1.1);
            xlim(ax, xr);
            in_band = f >= xr(1) & f <= xr(2);          % scale panel to its own
            pk = max(y(in_band));                         % in-band peak
            if isempty(pk) || ~(pk > 0), pk = 1; end
            ylim(ax, [0 pk*1.1]);
            set(ax, 'FontSize', 7); box(ax,'off');
            set(ax, 'XTick', [10 20 30]);   % x tick labels kept on every panel
            if rr == last_row(cc)
                xlabel(ax, 'Frequency (Hz)', 'FontSize', 8);
            end
        end
        if rr == 1, title(ax, pair_names{cc}, 'FontSize', 10); end
        if cc == 1
            % participant row label, drawn as a text child so it survives on
            % the empty (spine-only) Brain-EMG panels, where axis off would
            % otherwise hide a ylabel
            text(ax, -0.30, 0.5, sprintf('P%d', rr), 'Units','normalized', ...
                'Rotation', 90, 'HorizontalAlignment','center', ...
                'FontWeight','bold', 'FontSize', 10);
        end
    end
end
% one y-axis units label for the whole grid (units are x1e-3 throughout)
hax_lbl = axes('Parent', hfig_grid, 'Position', [0 0 1 1], 'Visible', 'off');
text(hax_lbl, 0.03, 0.5, 'Coherence (x1e-3)', 'Rotation', 90, ...
    'HorizontalAlignment', 'center', 'FontSize', 11);
if saveFigs
    grid_base = fullfile(fig_dir,['coherence_spectra_grid_raw' cfg.fig_suffix]);
    savefig(hfig_grid, [grid_base '.fig']);
    saveas(hfig_grid,  [grid_base '.png']);
    exportgraphics(hfig_grid, [grid_base '.pdf'], 'ContentType','vector');
    print(hfig_grid,   [grid_base '.svg'], '-dsvg', '-painters');
end

%% Peak coherence frequency summary (10-35 Hz) — one horizontal box per pair
% Brain pairs use the 7 common subjects; Cord-EMG uses all 9. Each participant
% is a labeled point; P-index is the participant's position in the canonical
% grid_subs order (brain subjects P1-P7, then spine-only P8-P9).
pf_BE = global_peak_freq(:,1);   % 7 x 1 (subs order)
pf_BC = global_peak_freq(:,2);
pf_SE = global_peak_freq(:,3);

pf_SE_all = nan(nRows,1);        % Cord-EMG peak freq, canonical order (n=9)
for cc = 1:nRows
    sub_c = grid_subs{cc};
    si = find(strcmp(subs, sub_c));
    if ~isempty(si)
        pf_SE_all(cc) = pf_SE(si);
    else
        ei = find(strcmp({spineEMG_extra.sub}, sub_c));
        if ~isempty(ei), pf_SE_all(cc) = spineEMG_extra(ei).peakf; end
    end
end

pf_sets  = {pf_BE, pf_BC, pf_SE_all};
pf_names = {'Brain-EMG','Brain-Cord','Cord-EMG'};
pf_cols  = {col_BE, col_BS, col_SE};

% per-participant peak-frequency matrix (spine P-index x pair), used to trace
% each participant across the 3 pairs; brain pairs are NaN for spine-only subs
pf_by_part = nan(nRows, 3);
for cc = 1:nRows
    sub_c = grid_subs{cc};
    si = find(strcmp(subs, sub_c));
    if ~isempty(si)
        pf_by_part(cc,1) = pf_BE(si);
        pf_by_part(cc,2) = pf_BC(si);
    end
    pf_by_part(cc,3) = pf_SE_all(cc);
end
% per-participant colours consistent with plot_freq_correlation_plots.m: the 7
% brain-cap-spine subjects take lines(nSubs) in the SAME subject order used
% there (subs = intersect(subs_brain,subs_spine,'stable')), so a given subject
% keeps its colour across scripts. The 2 spine-only subjects (not present in
% that script) get distinct extra colours that cannot collide with lines(7).
base_cmap  = lines(nSubs);
extra_cols = [0.00 0.00 0.00;    % OP00220 — near-black
              0.60 0.40 0.20];   % OP00226 — brown
pt_cmap = zeros(nRows,3);
ie = 0;
for cc = 1:nRows
    si = find(strcmp(subs, grid_subs{cc}));
    if ~isempty(si)
        pt_cmap(cc,:) = base_cmap(si,:);
    else
        ie = ie + 1;
        pt_cmap(cc,:) = extra_cols(min(ie,size(extra_cols,1)),:);
    end
end

hfig_pf = figure('Color','w','Position',[100 100 640 420]); hold on;
box_hh = 0.28;
for rr = 1:3
    v = pf_sets{rr}; v = v(isfinite(v));
    if isempty(v), continue; end
    q1 = quantile(v,0.25); q3 = quantile(v,0.75); md = median(v);
    wlo = min(v); whi = max(v); yc = rr;
    plot([wlo q1],[yc yc],'-','Color',[0.2 0.2 0.2],'LineWidth',1.2);
    plot([q3 whi],[yc yc],'-','Color',[0.2 0.2 0.2],'LineWidth',1.2);
    plot([wlo wlo],yc+[-0.1 0.1],'-','Color',[0.2 0.2 0.2],'LineWidth',1.2);
    plot([whi whi],yc+[-0.1 0.1],'-','Color',[0.2 0.2 0.2],'LineWidth',1.2);
    fill([q1 q3 q3 q1], yc+[-1 -1 1 1]*box_hh, [0.85 0.85 0.85], ...
        'FaceAlpha',0.35,'EdgeColor',[0.2 0.2 0.2],'LineWidth',1.2);
    plot([md md], yc+[-1 1]*box_hh,'-','Color',[0.2 0.2 0.2],'LineWidth',2.2);
end
% per-participant traces + markers, one colour per participant
rng(0);
pf_yj = 0.12*randn(nRows,1);   % one y-jitter per participant, shared across pairs
% P8 (OP00220, spine-only, lone Cord-EMG dot) sits exactly on P5 (OP00225) at
% ~19 Hz — nudge it up so the two markers/labels separate
pf_yj(strcmp(grid_subs,'OP00220')) = pf_yj(strcmp(grid_subs,'OP00220')) + 0.30;
for cc = 1:nRows
    yy = (1:3) + pf_yj(cc);
    xx = pf_by_part(cc,:);
    ok = isfinite(xx);
    if sum(ok) >= 2   % connecting trace needs at least two pairs
        plot(xx(ok), yy(ok), '-', 'Color', [pt_cmap(cc,:) 0.55], ...
            'LineWidth', 1.3, 'HandleVisibility','off');
    end
    for rr = 1:3
        if ~ok(rr), continue; end
        plot(xx(rr), yy(rr), 'o','MarkerSize',7, ...
            'MarkerFaceColor',pt_cmap(cc,:),'MarkerEdgeColor','w','LineWidth',0.8);
    end
    % label each participant at their Cord-EMG marker (no leader line)
    text(xx(3)+0.25, yy(3), sprintf('P%d', cc), 'FontSize',8, ...
        'Color', pt_cmap(cc,:)*0.75, 'FontWeight','bold','VerticalAlignment','middle');
end
xlim([fband(1)-1 fband(2)+1]); ylim([0.4 3.9]);
set(gca,'YTick',1:3,'YTickLabel',pf_names,'FontSize',12,'TickDir','out');
xlabel('Peak coherence frequency (Hz)','FontSize',12);
box off;
if saveFigs
    pf_base = fullfile(fig_dir,['peak_coherence_frequency_summary' cfg.fig_suffix]);
    savefig(hfig_pf, [pf_base '.fig']);
    saveas(hfig_pf,  [pf_base '.png']);
    exportgraphics(hfig_pf, [pf_base '.pdf'], 'ContentType','vector');
    print(hfig_pf,   [pf_base '.svg'], '-dsvg', '-painters');
end

fprintf('\n=== Peak coherence frequency (%.0f-%.0f Hz) ===\n', fband(1), fband(2));
fprintf('  %-12s  %9s  %9s  %9s\n','Sub','Brain-EMG','Brain-Cord','Cord-EMG');
for cc = 1:nRows
    sub_c = grid_subs{cc};
    si = find(strcmp(subs, sub_c));
    be = NaN; bc = NaN;
    if ~isempty(si), be = pf_BE(si); bc = pf_BC(si); end
    fprintf('  %-12s  %9.2f  %9.2f  %9.2f\n', sub_c, be, bc, pf_SE_all(cc));
end
fprintf('  %-12s  %9.2f  %9.2f  %9.2f\n','Median', ...
    median(pf_BE,'omitnan'), median(pf_BC,'omitnan'), median(pf_SE_all,'omitnan'));
fprintf('  %-12s  %9.2f  %9.2f  %9.2f\n','MAD', ...
    mad(pf_BE,1), mad(pf_BC,1), mad(pf_SE_all,1));
fprintf('  ---------------------------------------------\n');
fprintf('  Median peak frequency by pair (Hz):\n');
fprintf('    Brain-EMG  (n=%d): %.2f\n', sum(isfinite(pf_BE)),     median(pf_BE,'omitnan'));
fprintf('    Brain-Cord (n=%d): %.2f\n', sum(isfinite(pf_BC)),     median(pf_BC,'omitnan'));
fprintf('    Cord-EMG   (n=%d): %.2f\n', sum(isfinite(pf_SE_all)), median(pf_SE_all,'omitnan'));
fprintf('===============================================\n');

%% Brain-spine peak coherence vs threshold — boxplot
% Full coherence + NeuroSpec analytic CL. Box = median/IQR, whiskers to
% min/max, individual subjects jittered and colored by significance,
% P1 highlighted with an open ring.
conf95_BS_all = arrayfun(@(s) s.conf95_BS, results)';
ratio_BS      = peak_coh_brainSpine_full ./ conf95_BS_all;
pass_flag_BS  = ratio_BS > 1;
%% Save brain-spine boxplot data for figure regeneration
boxplot_data.peak_coh_brainSpine_full = peak_coh_brainSpine_full;
boxplot_data.conf95_BS_all            = conf95_BS_all;
boxplot_data.ratio_BS                 = ratio_BS;
boxplot_data.pass_flag_BS             = pass_flag_BS;
boxplot_data.subs                     = subs;
boxplot_data.p1_idx                   = p1_idx;
boxplot_data.nSubs                    = nSubs;
save(fullfile(cfg.save_dir, ...
    ['brainspine_boxplot_data' cfg.fig_suffix '.mat']), ...
    'boxplot_data');
fprintf('  Brain-spine boxplot data saved.\n');

col_pass = col_BS;
col_fail = [0.65 0.65 0.65];
col_edge = [0.15 0.15 0.15];
col_box  = [0.93 0.93 0.93];
col_line = [0.20 0.20 0.20];
col_thr  = [0.10 0.10 0.10];

q1  = quantile(ratio_BS, 0.25);
q3  = quantile(ratio_BS, 0.75);
med = median(ratio_BS);
wlo = min(ratio_BS);
whi = max(ratio_BS);

box_hw = 0.22;
cap_hw = 0.10;
jitter_BS = linspace(-0.16, 0.16, nSubs);

y_top = max([ratio_BS(:); 1]) * 1.15;
y_top = ceil(y_top*2)/2;

hfig_box = figure('Color','w','Position',[100 100 600 480]);
ax = axes('Parent', hfig_box); hold(ax,'on');

plot(ax, [0.3 1.7], [1 1], '--', 'Color', col_thr, 'LineWidth', 1.5, 'HandleVisibility','off');

plot(ax, [1 1], [q3 whi], '-', 'Color', col_line, 'LineWidth', 1.2, 'HandleVisibility','off');
plot(ax, [1 1], [wlo q1], '-', 'Color', col_line, 'LineWidth', 1.2, 'HandleVisibility','off');
plot(ax, [1-cap_hw 1+cap_hw], [whi whi], '-', 'Color', col_line, 'LineWidth', 1.2, 'HandleVisibility','off');
plot(ax, [1-cap_hw 1+cap_hw], [wlo wlo], '-', 'Color', col_line, 'LineWidth', 1.2, 'HandleVisibility','off');

fill(ax, [1-box_hw 1+box_hw 1+box_hw 1-box_hw], [q1 q1 q3 q3], col_box, ...
    'EdgeColor', col_line, 'LineWidth', 1.2, 'HandleVisibility','off');

plot(ax, [1-box_hw 1+box_hw], [med med], '-', 'Color', col_line, 'LineWidth', 2.2, 'HandleVisibility','off');

for ss = 1:nSubs
    if pass_flag_BS(ss), c = col_pass; else, c = col_fail; end
    plot(ax, 1+jitter_BS(ss), ratio_BS(ss), 'o', 'MarkerFaceColor', c, 'MarkerEdgeColor', col_edge, ...
        'MarkerSize', 9, 'LineWidth', 1, 'HandleVisibility','off');
end

plot(ax, 1+jitter_BS(p1_idx), ratio_BS(p1_idx), 'o', 'MarkerFaceColor', 'none', ...
    'MarkerEdgeColor', 'k', 'MarkerSize', 18, 'LineWidth', 2, 'HandleVisibility','off');

h_pass = plot(ax, nan, nan, 'o', 'MarkerFaceColor', col_pass, 'MarkerEdgeColor', col_edge, 'MarkerSize', 9, 'LineWidth', 1);
h_fail = plot(ax, nan, nan, 'o', 'MarkerFaceColor', col_fail, 'MarkerEdgeColor', col_edge, 'MarkerSize', 9, 'LineWidth', 1);
h_thr  = plot(ax, [NaN NaN], [NaN NaN], '--', 'Color', col_thr, 'LineWidth', 1.5);
legend(ax, [h_pass h_fail h_thr], {'Significant','Nonsignificant','Threshold'}, ...
    'Location','best', 'FontSize', 10, 'Box','off');

xlim(ax, [0.3 1.7]); ylim(ax, [0 y_top]);
set(ax, 'XTick', 1, 'XTickLabel', sprintf('Participants'));
set(ax, 'YTick', 0:0.5:y_top);
set(ax, 'FontSize', 11, 'FontName', 'Helvetica', 'TickDir', 'out', 'LineWidth', 1);
box(ax, 'off');

ylabel(ax, 'Peak coherence / threshold', 'FontSize', 13, 'FontWeight', 'bold');
title(ax, 'Brain-spinal cord coherence', 'FontSize', 15, 'FontWeight', 'bold', 'Color', [0.15 0.15 0.15]);

if saveFigs
    savefig(hfig_box, fullfile(fig_dir, ['brainspine_peak_vs_threshold_boxplot' cfg.fig_suffix '.fig']));
    saveas(hfig_box,  fullfile(fig_dir, ['brainspine_peak_vs_threshold_boxplot' cfg.fig_suffix '.png']));
end

%% Power statistics
fprintf('\n=== Group: paired t-tests %.0f-%.0f Hz band power ===\n', fband(1), fband(2));
okB = isfinite(Pstat_brain) & isfinite(Prest_brain);
[~,pb,~,sb] = ttest(log(Pstat_brain(okB)), log(Prest_brain(okB)));
fprintf('  Brain: n=%d, t(%d)=%.3f, p=%.4g\n', sum(okB), sb.df, sb.tstat, pb);
okS = isfinite(Pstat_spine) & isfinite(Prest_spine);
[~,ps,~,ss_] = ttest(log(Pstat_spine(okS)), log(Prest_spine(okS)));
fprintf('  Spine: n=%d, t(%d)=%.3f, p=%.4g\n', sum(okS), ss_.df, ss_.tstat, ps);

fprintf('\n=== Participant 1 — trial-level paired t-test (%.0f-%.0f Hz power) ===\n', fband(1), fband(2));
if ~isempty(p1_brain_stat_trials)
    okB1 = isfinite(p1_brain_stat_trials) & isfinite(p1_brain_rest_trials) & ...
           p1_brain_stat_trials > 0 & p1_brain_rest_trials > 0;
    [~,pb1,~,sb1] = ttest(log(p1_brain_stat_trials(okB1)), log(p1_brain_rest_trials(okB1)));
    fprintf('  Brain VE: n=%d trials, t(%d)=%.3f, p=%.4g\n', sum(okB1), sb1.df, sb1.tstat, pb1);
    okS1 = isfinite(p1_spine_stat_trials) & isfinite(p1_spine_rest_trials) & ...
           p1_spine_stat_trials > 0 & p1_spine_rest_trials > 0;
    [~,ps1,~,ss1] = ttest(log(p1_spine_stat_trials(okS1)), log(p1_spine_rest_trials(okS1)));
    fprintf('  Spine VE: n=%d trials, t(%d)=%.3f, p=%.4g\n', sum(okS1), ss1.df, ss1.tstat, ps1);
else
    fprintf('  P1 trial data not available.\n');
end
fprintf('=======================================================\n');

%% FOOOF
fprintf('\n=== Group: FOOOF periodic power (%.0f-%.0f Hz) ===\n', fband(1), fband(2));
fprintf('  Sub           Brain_stat   Brain_rest   Spine_stat   Spine_rest\n');
for ss = 1:nSubs
    fprintf('  %-12s  %10.4e   %10.4e   %10.4e   %10.4e\n', subs{ss}, ...
        Pstat_brain_fooof(ss), Prest_brain_fooof(ss), ...
        Pstat_spine_fooof(ss), Prest_spine_fooof(ss));
end
fprintf('  Median stat:  %10.4e              %10.4e\n', ...
    median(Pstat_brain_fooof,'omitnan'), median(Pstat_spine_fooof,'omitnan'));
fprintf('  Median rest:  %10.4e              %10.4e\n', ...
    median(Prest_brain_fooof,'omitnan'), median(Prest_spine_fooof,'omitnan'));
okBf = Pstat_brain_fooof>0 & Prest_brain_fooof>0 & isfinite(Pstat_brain_fooof) & isfinite(Prest_brain_fooof);
okSf = Pstat_spine_fooof>0 & Prest_spine_fooof>0 & isfinite(Pstat_spine_fooof) & isfinite(Prest_spine_fooof);
fprintf('\n  Paired t-test (log, zeros excluded):\n');
if sum(okBf)>=2
    [~,pbf,~,sbf] = ttest(log(Pstat_brain_fooof(okBf)), log(Prest_brain_fooof(okBf)));
    fprintf('  Brain: n=%d/%d, t(%d)=%.3f, p=%.4g\n', sum(okBf), nSubs, sbf.df, sbf.tstat, pbf);
else
    fprintf('  Brain: insufficient pairs (n=%d)\n', sum(okBf));
end
if sum(okSf)>=2
    [~,psf,~,ssf] = ttest(log(Pstat_spine_fooof(okSf)), log(Prest_spine_fooof(okSf)));
    fprintf('  Spine: n=%d/%d, t(%d)=%.3f, p=%.4g\n', sum(okSf), nSubs, ssf.df, ssf.tstat, psf);
else
    fprintf('  Spine: insufficient pairs (n=%d)\n', sum(okSf));
end
fprintf('=======================================================\n');

%% P1 beta band power
fprintf('\n=== Participant 1 — beta band power (%.0f-%.0f Hz) ===\n', fband(1), fband(2));
fprintf('  Brain VE:    static=%.4e  rest=%.4e  ratio=%.3f\n', ...
    Pstat_brain(p1_idx), Prest_brain(p1_idx), Pstat_brain(p1_idx)/Prest_brain(p1_idx));
fprintf('  Spine VE:    static=%.4e  rest=%.4e  ratio=%.3f\n', ...
    Pstat_spine(p1_idx), Prest_spine(p1_idx), Pstat_spine(p1_idx)/Prest_spine(p1_idx));
fprintf('  Brain FOOOF: static=%.4e  rest=%.4e\n', ...
    Pstat_brain_fooof(p1_idx), Prest_brain_fooof(p1_idx));
fprintf('  Spine FOOOF: static=%.4e  rest=%.4e\n', ...
    Pstat_spine_fooof(p1_idx), Prest_spine_fooof(p1_idx));
fprintf('====================================================\n');

%% Supplementary table: peak coherence magnitude
tbl_vars = {'BrainEMG_MT','BrainSpine_MT','SpineEMG_MT', ...
            'BrainEMG_NT','BrainSpine_NT','SpineEMG_NT'};
% Row labels as P#/subject (canonical numbering) so readers can map the
% participant labels used in the figures back to subject IDs.
row_lbl_main = arrayfun(@(k) sprintf('P%d/%s', k, subs{k}), ...
    1:numel(subs), 'UniformOutput', false);
peak_coh_table_full = array2table(...
    [peak_coh_brainEMG, peak_coh_brainSpine, peak_coh_spineEMG, ...
     peak_coh_brainEMG_nt, peak_coh_brainSpine_nt, peak_coh_spineEMG_nt], ...
    'VariableNames',tbl_vars,'RowNames',row_lbl_main);

% Append the cord-EMG-only subjects (P8/P9): brain pairs are NaN, SpineEMG
% filled. They follow the 7 brain subjects, matching the canonical numbering.
if ~isempty(spineEMG_extra)
    ex_mt   = [spineEMG_extra.peak_mt]';
    ex_nt   = [spineEMG_extra.peak_nt]';
    nEx     = numel(ex_mt);
    row_lbl_ex = arrayfun(@(j) sprintf('P%d/%s', numel(subs)+j, spineEMG_extra(j).sub), ...
        1:nEx, 'UniformOutput', false);
    ex_rows = array2table([nan(nEx,2), ex_mt, nan(nEx,2), ex_nt], ...
        'VariableNames',tbl_vars,'RowNames',row_lbl_ex);
    peak_coh_table_full = [peak_coh_table_full; ex_rows];
end

% SpineEMG (Cord-EMG) summary spans all available subjects (n=9); the two
% brain pairs stay n=7. mad() has no omitnan, so filter finite values.
spineEMG_mt_all = peak_coh_table_full.SpineEMG_MT;
spineEMG_nt_all = peak_coh_table_full.SpineEMG_NT;
madfin = @(x) mad(x(isfinite(x)),1);
coh_summary = [median(peak_coh_brainEMG,'omitnan'),    median(peak_coh_brainSpine,'omitnan'),    median(spineEMG_mt_all,'omitnan'), ...
    median(peak_coh_brainEMG_nt,'omitnan'), median(peak_coh_brainSpine_nt,'omitnan'), median(spineEMG_nt_all,'omitnan'); ...
    mad(peak_coh_brainEMG,1),    mad(peak_coh_brainSpine,1),    madfin(spineEMG_mt_all), ...
    mad(peak_coh_brainEMG_nt,1), mad(peak_coh_brainSpine_nt,1), madfin(spineEMG_nt_all)];
supp_table_1b = [peak_coh_table_full; array2table(coh_summary, ...
    'VariableNames',tbl_vars,'RowNames',{'Median','MAD'})];
fprintf('\n  Supplementary Table 1b: Peak coherence (R2)\n'); disp(supp_table_1b);
try
    writetable(supp_table_1b, fullfile(fig_dir,'SuppTable1b_peak_coherence_BS.csv'),'WriteRowNames',true);
catch ME, warning('Could not save Table 1b: %s', ME.message); end

%% Collect PSI
for ss = 1:nSubs
    if isfield(results(ss),'psi_brainEMG')
        psi_brainEMG(ss)   = results(ss).psi_brainEMG;
        psi_brainSpine(ss) = results(ss).psi_brainSpine;
        psi_spineEMG(ss)   = results(ss).psi_spineEMG;
    end
end

%% Directionality comparison (Halliday R2 vs PSI)
pair_labels_dir = {'Brain-EMG','Brain-Spine','Spine-EMG'};
hal_fwd = nan(nSubs,3); hal_rev = nan(nSubs,3);
for ss = 1:nSubs
    hal_fwd(ss,:) = [results(ss).brainEMG.forward_area, results(ss).brainSpine.forward_area, results(ss).spineEMG.forward_area];
    hal_rev(ss,:) = [results(ss).brainEMG.reverse_area, results(ss).brainSpine.reverse_area, results(ss).spineEMG.reverse_area];
end
psi_mat = [psi_brainEMG, psi_brainSpine, psi_spineEMG];
hfig_dir_comp = figure('Color','w','Position',[100 100 900 560]);
for pp = 1:3
    jit = 0.08;
    subplot(2,3,pp); hold on;
    for ss = 1:nSubs
        xj = [1+jit*randn, 2+jit*randn]; yj = [hal_fwd(ss,pp), hal_rev(ss,pp)];
        if ss==p1_idx
            plot(xj,yj,'-o','Color',[0.1 0.1 0.1],'LineWidth',2.5,'MarkerSize',8, ...
                'MarkerFaceColor',[0.1 0.1 0.1],'MarkerEdgeColor','w');
        else
            plot(xj,yj,'-','Color',[0.75 0.75 0.75],'LineWidth',0.8);
            scatter(xj(1),yj(1),35,[0.2 0.4 0.8],'filled','MarkerFaceAlpha',0.55,'MarkerEdgeColor','none');
            scatter(xj(2),yj(2),35,[0.9 0.3 0.3],'filled','MarkerFaceAlpha',0.55,'MarkerEdgeColor','none');
        end
    end
    plot(1,median(hal_fwd(:,pp),'omitnan'),'s','MarkerSize',12,'MarkerFaceColor',[0.2 0.4 0.8],'MarkerEdgeColor','k','LineWidth',1.5);
    plot(2,median(hal_rev(:,pp),'omitnan'),'s','MarkerSize',12,'MarkerFaceColor',[0.9 0.3 0.3],'MarkerEdgeColor','k','LineWidth',1.5);
    set(gca,'XTick',[1 2],'XTickLabel',{'Forward','Reverse'},'FontSize',14);
    ylabel('Coherence area','FontSize',14);
    title(sprintf('%s\nHalliday R2',pair_labels_dir{pp}),'FontSize',14);
    ylim([0 max([hal_fwd(:,pp);hal_rev(:,pp)])*1.3]); xlim([0.5 2.5]); box off;
    subplot(2,3,3+pp); hold on;
    yline(0,'-k','LineWidth',1,'HandleVisibility','off');
    for ss = 1:nSubs
        xj = 1+jit*randn; yj = psi_mat(ss,pp);
        if isnan(yj), continue; end
        if ss==p1_idx
            plot(xj,yj,'o','Color',[0.1 0.1 0.1],'LineWidth',2.5,'MarkerSize',10, ...
                'MarkerFaceColor',[0.1 0.1 0.1],'MarkerEdgeColor','w');
        else
            scatter(xj,yj,45,[0.2 0.7 0.3],'filled','MarkerFaceAlpha',0.55,'MarkerEdgeColor','none');
        end
    end
    med_psi = median(psi_mat(:,pp),'omitnan'); mad_psi = mad(psi_mat(:,pp),1);
    plot(1,med_psi,'s','MarkerSize',12,'MarkerFaceColor',[0.2 0.7 0.3],'MarkerEdgeColor','k','LineWidth',1.5);
    errorbar(1,med_psi,mad_psi,'k','LineWidth',1.5,'CapSize',8);
    set(gca,'XTick',1,'XTickLabel',{'PSI'},'FontSize',14);
    ylabel('Summed PSI (beta band)','FontSize',14);
    title(sprintf('%s\nPSI',pair_labels_dir{pp}),'FontSize',14);
    psi_lim = max(abs(psi_mat(:,pp)),[],'omitnan')*1.4;
    if isnan(psi_lim)||psi_lim==0, psi_lim=1; end
    ylim([-psi_lim psi_lim]); xlim([0.5 1.5]); box off;
end
sgtitle('Directionality: Halliday R2 (top) vs PSI (bottom, BS)','FontSize',14,'Interpreter','none');
if saveFigs
    savefig(hfig_dir_comp, fullfile(fig_dir,['group_directionality_comparison' cfg.fig_suffix '.fig']));
    saveas(hfig_dir_comp,  fullfile(fig_dir,['group_directionality_comparison' cfg.fig_suffix '.png']));
end

%% Peak latency plot + printout
brainEMG_lat   = abs(arrayfun(@(s) s.brainEMG.peak_latency,   results));
brainSpine_lat = abs(arrayfun(@(s) s.brainSpine.peak_latency, results));
spineEMG_lat   = abs(arrayfun(@(s) s.spineEMG.peak_latency,   results));
lightGrey=[0.75 0.75 0.75]; darkGrey=[0.25 0.25 0.25]; medianGrey=[0.10 0.10 0.10];
x_lat=[1 2 3];
hfig_lat = figure('Color','w'); hold on; hP1=[]; hMed=[];
for ss = 1:nSubs
    y = [brainEMG_lat(ss), brainSpine_lat(ss), spineEMG_lat(ss)];
    if any(isnan(y)), continue; end
    if ss==p1_idx
        hP1 = plot(x_lat,y,'-o','LineWidth',3.5,'MarkerSize',9,'Color',darkGrey, ...
            'MarkerFaceColor',darkGrey,'MarkerEdgeColor','w','DisplayName',subs{ss});
    else
        plot(x_lat,y,'-o','Color',lightGrey,'LineWidth',0.8,'MarkerSize',6, ...
            'MarkerFaceColor',lightGrey,'MarkerEdgeColor','w');
    end
end
scatter(ones(nSubs,1),   brainEMG_lat',   60,lightGrey,'filled','MarkerFaceAlpha',0.5,'MarkerEdgeColor','none');
scatter(2*ones(nSubs,1), brainSpine_lat', 60,lightGrey,'filled','MarkerFaceAlpha',0.5,'MarkerEdgeColor','none');
scatter(3*ones(nSubs,1), spineEMG_lat',   60,lightGrey,'filled','MarkerFaceAlpha',0.5,'MarkerEdgeColor','none');
med_vals=[median(brainEMG_lat,'omitnan'),median(brainSpine_lat,'omitnan'),median(spineEMG_lat,'omitnan')];
mad_vals=[mad(brainEMG_lat,1),mad(brainSpine_lat,1),mad(spineEMG_lat,1)];
hMed=plot(x_lat,med_vals,'--s','LineWidth',2.5,'MarkerSize',9,'Color',medianGrey, ...
    'MarkerFaceColor','w','MarkerEdgeColor',medianGrey);
errorbar(x_lat,med_vals,mad_vals,'Color',medianGrey,'LineStyle','none','LineWidth',1.5,'CapSize',10);
leg_handles=[]; leg_labels={};
if ~isempty(hP1), leg_handles=[leg_handles,hP1]; leg_labels=[leg_labels,{'Participant 1'}]; end
leg_handles=[leg_handles,hMed]; leg_labels=[leg_labels,{'Median +/- MAD'}];
legend(leg_handles,leg_labels,'Location','northwest','Box','off','FontSize',14);
xlim([0.5 3.8]); xticks(x_lat);
xticklabels({'Brain<->EMG','Brain<->Spine','Spine<->EMG'});
ylabel('Peak latency (ms)','FontSize',14); set(gca,'FontSize',14); grid on; box on;
title('Cross-correlation peak latencies (BS)','Interpreter','none','FontSize',14);
if saveFigs
    savefig(hfig_lat, fullfile(fig_dir,['group_peak_latencies' cfg.fig_suffix '.fig']));
    saveas(hfig_lat,  fullfile(fig_dir,['group_peak_latencies' cfg.fig_suffix '.png']));
end

fprintf('\n=== Peak latencies (ms) ===\n');
fprintf('  %-20s  %8s  %8s  %8s\n', '', 'Brain-EMG', 'Br-Spine', 'Spine-EMG');
fprintf('  %-20s  %8.1f  %8.1f  %8.1f\n', 'P1', ...
    brainEMG_lat(p1_idx), brainSpine_lat(p1_idx), spineEMG_lat(p1_idx));
fprintf('  %-20s  %8.1f  %8.1f  %8.1f\n', 'Group median', med_vals(1), med_vals(2), med_vals(3));
fprintf('  %-20s  %8.1f  %8.1f  %8.1f\n', 'Group MAD',    mad_vals(1), mad_vals(2), mad_vals(3));
fprintf('  Individual subjects:\n');
for ss = 1:nSubs
    fprintf('  %-20s  %8.1f  %8.1f  %8.1f\n', subs{ss}, ...
        brainEMG_lat(ss), brainSpine_lat(ss), spineEMG_lat(ss));
end
fprintf('=======================================================\n');

%% Final summary: median frequency of peak coherence, per signal pair
fprintf('\n=== Median frequency of peak coherence by pair (%.0f-%.0f Hz) ===\n', fband(1), fband(2));
fprintf('  Brain-EMG  (n=%d): %.2f Hz\n', sum(isfinite(pf_BE)),     median(pf_BE,'omitnan'));
fprintf('  Brain-Cord (n=%d): %.2f Hz\n', sum(isfinite(pf_BC)),     median(pf_BC,'omitnan'));
fprintf('  Cord-EMG   (n=%d): %.2f Hz\n', sum(isfinite(pf_SE_all)), median(pf_SE_all,'omitnan'));
fprintf('===============================================================\n');

end  % run_coherence_spectra


%% =========================================================================
%  UTILITY FUNCTIONS
%% =========================================================================

function fname = get_datafile(data_root, sub)
run = '001'; if strcmp(sub,'OP00224'), run = '002'; end
fname = fullfile(data_root, ['sub-' sub], 'ses-001', 'meg', ...
    sprintf('pmergedoe1000mspddfflo45hi45hfcstatic_%s_array1.mat', run));
end

function lats = get_top_peaks(t, lat_min, lat_max, nPeaks)
valid_idx = t(:,1) >= lat_min & t(:,1) <= lat_max;
rho_vals  = t(valid_idx, 3); lag_vals = t(valid_idx, 1);
[~, locs] = findpeaks(abs(rho_vals));
if isempty(locs), locs = (1:numel(rho_vals))'; end
locs = locs(abs(lag_vals(locs)) > 2);
if isempty(locs), [~,locs] = max(abs(rho_vals)); end
[~,ord] = sort(abs(rho_vals(locs)),'descend');
locs_sorted = locs(ord);
lats = lag_vals(locs_sorted(1:min(nPeaks,numel(locs_sorted))));
end

function [peak_rho, peak_latency] = select_peak(t, lat_min, lat_max, zero_excl_ms)
valid_idx = t(:,1) >= lat_min & t(:,1) <= lat_max;
rho_vals  = t(valid_idx,3); lag_vals = t(valid_idx,1);
if nargin < 4 || isempty(zero_excl_ms), zero_excl_ms = 2; end
[~, locs] = findpeaks(abs(rho_vals));
if isempty(locs), locs = (1:numel(rho_vals))'; end
locs = locs(abs(lag_vals(locs)) > zero_excl_ms);
if isempty(locs)
    [~,idx_max] = max(abs(rho_vals));
    peak_rho = rho_vals(idx_max); peak_latency = lag_vals(idx_max);
    fprintf('All peaks within +/-%g ms. Returning global max.\n', zero_excl_ms); return
end
[~,ord] = sort(abs(rho_vals(locs)),'descend');
locs_sorted = locs(ord);
nPrint = min(5, numel(locs_sorted));
fprintf('\nTop %d peaks outside +/-%g ms:\n', nPrint, zero_excl_ms);
fprintf('  Rank   Lag (ms)     rho        |rho|\n');
for k = 1:nPrint
    ii = locs_sorted(k);
    fprintf('   %d     %7.2f    %+8.4g    %8.4g\n', k, lag_vals(ii), rho_vals(ii), abs(rho_vals(ii)));
end
idx_max = locs_sorted(1);
peak_rho = rho_vals(idx_max); peak_latency = lag_vals(idx_max);
end

function [partial_coh, freq_vec, full_coh] = mt_partial_coherence(x, y, z, samp_rate, seg_len, NW)
K = round(2*NW) - 1;
[tapers, ~] = dpss(seg_len, NW, K);
n_segs  = floor(length(x) / seg_len);
n_freqs = seg_len/2 + 1;
f_xx=zeros(n_freqs,1); f_yy=zeros(n_freqs,1); f_zz=zeros(n_freqs,1);
f_xy=zeros(n_freqs,1); f_xz=zeros(n_freqs,1); f_yz=zeros(n_freqs,1);
for seg = 1:n_segs
    idx = (seg-1)*seg_len + (1:seg_len);
    xs=x(idx)-mean(x(idx)); ys=y(idx)-mean(y(idx)); zs=z(idx)-mean(z(idx));
    for k = 1:K
        tap=tapers(:,k);
        Xk=fft(xs.*tap); Xk=Xk(1:n_freqs);
        Yk=fft(ys.*tap); Yk=Yk(1:n_freqs);
        Zk=fft(zs.*tap); Zk=Zk(1:n_freqs);
        f_xx=f_xx+abs(Xk).^2; f_yy=f_yy+abs(Yk).^2; f_zz=f_zz+abs(Zk).^2;
        f_xy=f_xy+Xk.*conj(Yk); f_xz=f_xz+Xk.*conj(Zk); f_yz=f_yz+Yk.*conj(Zk);
    end
end
denom=n_segs*K;
f_xx=f_xx/denom; f_yy=f_yy/denom; f_zz=f_zz/denom;
f_xy=f_xy/denom; f_xz=f_xz/denom; f_yz=f_yz/denom;
f_xy_z = f_xy - (f_xz.*conj(f_yz))./f_zz;
f_xx_z = f_xx - (abs(f_xz).^2)./f_zz;
f_yy_z = f_yy - (abs(f_yz).^2)./f_zz;
partial_coh = abs(f_xy_z).^2./(f_xx_z.*f_yy_z);
partial_coh = max(0,min(1,real(partial_coh)));
full_coh    = abs(f_xy).^2./(f_xx.*f_yy);
full_coh    = max(0,min(1,real(full_coh)));
partial_coh(1)=0; full_coh(1)=0;
freq_vec = (0:n_freqs-1)'*(samp_rate/seg_len);
end

function s = pfmt_short(p)
if p < 0.001
    s = '<.001';
else
    s = sprintf('%.3f', p);
end
end

function s = conditional_str(cond, str_true, str_false)
if cond, s = str_true; else, s = str_false; end
end
