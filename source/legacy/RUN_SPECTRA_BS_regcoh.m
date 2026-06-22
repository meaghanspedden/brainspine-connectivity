%% RUN_SPECTRA_BS.m
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
    
    cfg.lag_mode = 'nonzero';   % 'nonzero' or 'full' — see header comment
    
    cfg.run_surrogates  = false;    % set false to skip (uses NeuroSpec CL as fallback)
    cfg.n_surrogates    = 500;
    cfg.surr_pctile_BE  = 95;
    cfg.surr_pctile_BS  = 95;     
    cfg.surr_pctile_SE  = 95;
    
    cfg.fband         = [10 35];   % analysis band (significance, power, coherence, PSI)
    
    cfg.seg_pwr = 11;
    cfg.lat_min = -50;
    cfg.lat_max =  50;
    
    saveFigs     = 1;
    cfg.saveFigs = saveFigs;
    
    assert(any(strcmp(cfg.lag_mode, {'full','nonzero'})), ...
        'cfg.lag_mode must be ''full'' or ''nonzero'', got ''%s''', cfg.lag_mode);
    
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
    fprintf('  Lag mode: %s\n', cfg.lag_mode)
    fprintf('  Analysis band:  %.0f-%.0f Hz\n', cfg.fband(1), cfg.fband(2))
    fprintf('  Run surrogates: %d\n', cfg.run_surrogates)
    fprintf('  Surr percentiles: BE=%d  BS=%d  SE=%d\n', ...
        cfg.surr_pctile_BE, cfg.surr_pctile_BS, cfg.surr_pctile_SE)
    fprintf('===========================\n\n')
    
    if ~cfg.run_surrogates && strcmp(cfg.lag_mode, 'nonzero')
        warning(['Using the analytic NeuroSpec confidence limit as the significance threshold, but that ' ...
            'limit (Rosenberg et al. 1989) is derived for the full (all-lag) coherence statistic, not the ' ...
            'forward+reverse non-zero-lag decomposition selected via lag_mode=''nonzero''. Treat thresholds ' ...
            'from this fallback as approximate; set cfg.run_surrogates = true for a properly-matched null.']);
    end
    
    fprintf('\n>>> STEP B: Coherence spectra\n')
    for vi = 1:size(ve_configs,1)
        cfg.spine_ve_pattern = ve_configs{vi,2};
        cfg.spine_ve_varname = ve_configs{vi,3};
        cfg.fig_suffix       = sprintf('_%s_%s', ve_configs{vi,1}, cfg.lag_mode);
        cfg.fig_dir = fullfile(cfg.save_dir, 'figures', sprintf('spectra_%s_%s', ve_configs{vi,1}, cfg.lag_mode));
        if ~exist(cfg.fig_dir,'dir'), mkdir(cfg.fig_dir); end
        fprintf('\n>>> Running VE config: %s (lag mode: %s)\n', ve_configs{vi,1}, cfg.lag_mode);
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
    lag_mode        = cfg.lag_mode;
    n_surrogates    = cfg.n_surrogates;
    run_surrogates  = cfg.run_surrogates;
    surr_pctile_BE  = cfg.surr_pctile_BE;
    surr_pctile_BS  = cfg.surr_pctile_BS;
    surr_pctile_SE  = cfg.surr_pctile_SE;
    
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
    
        %% Surrogate thresholds
        nF = size(f1_mt, 1);
        if run_surrogates
            fprintf('  Computing surrogate thresholds (%d surrogates)...\n', n_surrogates);
            results(ss).surr_thresh_BE = compute_surrogate_threshold(...
                statBcont, statEMGcont, samp_rate, seg_pwr, n_surrogates, opt_str, surr_pctile_BE, lag_mode);
            results(ss).surr_thresh_BS = compute_surrogate_threshold(...
                statBcont, statScont,   samp_rate, seg_pwr, n_surrogates, opt_str, surr_pctile_BS, lag_mode);
            results(ss).surr_thresh_SE = compute_surrogate_threshold(...
                statScont, statEMGcont, samp_rate, seg_pwr, n_surrogates, opt_str, surr_pctile_SE, lag_mode);
            freq_ax_tmp = f1_mt(:,1)';
            bmask_tmp   = freq_ax_tmp >= fband(1) & freq_ax_tmp <= fband(2);
            fprintf('  Surr thresholds (mean in band): BE=%.4f  BS=%.4f  SE=%.4f\n', ...
                mean(results(ss).surr_thresh_BE(bmask_tmp)), ...
                mean(results(ss).surr_thresh_BS(bmask_tmp)), ...
                mean(results(ss).surr_thresh_SE(bmask_tmp)));
        else
            results(ss).surr_thresh_BE = repmat(cl1_mt.ch_c95, 1, nF);
            results(ss).surr_thresh_BS = repmat(cl2_mt.ch_c95, 1, nF);
            results(ss).surr_thresh_SE = repmat(cl3_mt.ch_c95, 1, nF);
            fprintf('  Surrogates skipped — using NeuroSpec CL as threshold.\n');
        end
    
        %% P1 coherence spectra — standalone figure
        if strcmp(sub, 'OP00212')
            hfig_p1 = figure('Color','w','Position',[100 100 900 320]);
            fax_p1  = f1_mt(:,1);
            coh_BE  = get_coh(f1_mt, lag_mode)*1e3;
            coh_BS  = get_coh(f2_mt, lag_mode)*1e3;
            coh_SE  = get_coh(f3_mt, lag_mode)*1e3;
            pairs      = {coh_BE, coh_BS, coh_SE};
            cols       = {col_BE, col_BS, col_SE};
            titles_p1  = {'Brain-EMG','Brain-Spine','Spine-EMG'};
            surr_thresh_p1 = {results(1).surr_thresh_BE*1e3, ...
                              results(1).surr_thresh_BS*1e3, ...
                              results(1).surr_thresh_SE*1e3};
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
            sgtitle(sprintf('Participant 1 — coherence spectra (BSLaw, %s lag)', lag_mode),'FontSize',13,'Interpreter','none');
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
        coh_BE_full = get_coh(f1_mt, lag_mode);
        coh_BS_full = get_coh(f2_mt, lag_mode);
        coh_SE_full = get_coh(f3_mt, lag_mode);
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
    
        %% NT peak coherence (analysis band)
        fax_nt_pk = f1_nt(:,1)';
        band_nt   = fax_nt_pk >= fband(1) & fax_nt_pk <= fband(2);
        coh_BE_nt = get_coh(f1_nt, lag_mode);
        coh_BS_nt = get_coh(f2_nt, lag_mode);
        coh_SE_nt = get_coh(f3_nt, lag_mode);
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
        % NOT controlled by lag_mode.
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
    
    p1_idx = 1;  % Participant 1, highlighted throughout the group figures below
    
    %% =========================================================================
    %  FIGURES
    %% =========================================================================
    sig_BE   = nan(nSubs, 1);
    sig_BS   = nan(nSubs, 1);
    sig_SE   = nan(nSubs, 1);
    
    for ss = 1:nSubs
        fax = results(ss).freq_axis_mt;
        if isempty(fax), continue; end
        band_mask = fax >= fband(1) & fax <= fband(2);
        sig_BE(ss) = any(results(ss).coh_BE(band_mask) > results(ss).surr_thresh_BE(band_mask));
        sig_BS(ss) = any(results(ss).coh_BS(band_mask) > results(ss).surr_thresh_BS(band_mask));
        sig_SE(ss) = any(results(ss).coh_SE(band_mask) > results(ss).surr_thresh_SE(band_mask));
    end
    
    fprintf('\n=== Beta band coherence significance (surrogate threshold, lag mode: %s) ===\n', lag_mode);
    fprintf('  Sub           BE_sig  BS_sig  SE_sig\n');
    for ss = 1:nSubs
        fprintf('  %-12s  %s       %s       %s\n', subs{ss}, ...
            conditional_str(sig_BE(ss)==1,'YES','no '), ...
            conditional_str(sig_BS(ss)==1,'YES','no '), ...
            conditional_str(sig_SE(ss)==1,'YES','no '));
    end
    fprintf('  Significant BE: %d/%d\n', sum(sig_BE==1,'omitnan'), nSubs);
    fprintf('  Significant BS: %d/%d\n', sum(sig_BS==1,'omitnan'), nSubs);
    fprintf('  Significant SE: %d/%d\n', sum(sig_SE==1,'omitnan'), nSubs);
    fprintf('=============================================================\n');
    
    %% Partial vs full coherence — paired test across participants
    % Same comparison the aligned figure below shows visually (full vs
    % partial Brain-EMG coherence, controlling for spine), but reduced to
    % one paired value per subject: each subject's own full-coherence peak
    % (within fband) vs the partial coherence at that same frequency bin.
    % Both come from mt_partial_coherence (same estimator/tapering), so the
    % comparison isn't confounded by using two different coherence methods.
    % This is a separate technique from sp2a2_R2_mt and is unaffected by
    % lag_mode.
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

    %% Aligned partial coherence
    rel_win_p  = 12; rel_res_p = 0.25;
    rel_axis_p = (-rel_win_p : rel_res_p : rel_win_p);
    full_amat    = nan(nSubs, length(rel_axis_p));
    partial_amat = nan(nSubs, length(rel_axis_p));
    for ss = 1:nSubs
        fax    = results(ss).freq_axis_mt;
        full_v = results(ss).full_coh_mt_BE;
        part_v = results(ss).partial_mt_BE_spine;
        if isempty(fax), continue; end
        mask_full = fax >= fband(1) & fax <= fband(2);
        [~,ig] = max(full_v(mask_full));
        fc_s = fax(mask_full); fc_s = fc_s(ig);
        fax_rel  = fax - fc_s;
        in_range = rel_axis_p >= min(fax_rel) & rel_axis_p <= max(fax_rel);
        if sum(in_range) < 2, continue; end
        full_amat(ss,in_range)    = interp1(fax_rel, full_v, rel_axis_p(in_range), 'linear');
        partial_amat(ss,in_range) = interp1(fax_rel, part_v, rel_axis_p(in_range), 'linear');
    end
    mean_full    = mean(full_amat,1,'omitnan');
    mean_partial = mean(partial_amat,1,'omitnan');
    sem_full     = std(full_amat,0,1,'omitnan')/sqrt(nSubs);
    sem_partial  = std(partial_amat,0,1,'omitnan')/sqrt(nSubs);
    ok_cols      = ~all(isnan(full_amat),1);
    hfig_part_aligned = figure('Color','w','Position',[100 100 500 380]); hold on;
    fill([rel_axis_p(ok_cols) fliplr(rel_axis_p(ok_cols))], ...
        [(mean_full(ok_cols)+sem_full(ok_cols)) fliplr(mean_full(ok_cols)-sem_full(ok_cols))]*1e3, ...
        col_full,'FaceAlpha',0.2,'EdgeColor','none');
    fill([rel_axis_p(ok_cols) fliplr(rel_axis_p(ok_cols))], ...
        [(mean_partial(ok_cols)+sem_partial(ok_cols)) fliplr(mean_partial(ok_cols)-sem_partial(ok_cols))]*1e3, ...
        col_partial,'FaceAlpha',0.2,'EdgeColor','none');
    plot(rel_axis_p, mean_full*1e3,    '-', 'Color',col_full,    'LineWidth',2.5);
    plot(rel_axis_p, mean_partial*1e3, '--','Color',col_partial, 'LineWidth',2.5);
    xline(0,'--k','LineWidth',1.2,'HandleVisibility','off');
    xlim([-rel_win_p rel_win_p]);
    xlabel('Hz from Brain-EMG peak','FontSize',13); ylabel('Coherence (x1e-3)','FontSize',13);
    hf2 = plot(nan,nan,'-','Color',col_full,'LineWidth',2);
    hp2 = plot(nan,nan,'--','Color',col_partial,'LineWidth',2);
    legend([hf2 hp2],{'Full (BE)','Partial (spine removed)'},'Location','northeast','Box','off','FontSize',11);
    title('Brain-EMG: full vs partial coherence (aligned, BS)','FontSize',13);
    yl_pa = ylim;
    text(-rel_win_p*0.95, yl_pa(2)*0.92, ...
        sprintf('Wilcoxon p=%s  |  paired t p=%s  (n=%d)', pfmt_short(p_part_signrank), pfmt_short(p_part_ttest), n_part), ...
        'FontSize',10,'Color',[0.3 0.3 0.3]);
    set(gca,'FontSize',13); box off;
    if saveFigs
        savefig(hfig_part_aligned, fullfile(fig_dir,['partial_coherence_aligned' cfg.fig_suffix '.fig']));
    end
    
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
    sgtitle(sprintf('Group coherence spectra — aligned to each pair''s own peak (%s lag)', lag_mode),'FontSize',13);
    if saveFigs
        savefig(hfig_aligned, fullfile(fig_dir,['group_coherence_spectra_global' cfg.fig_suffix '.fig']));
        saveas(hfig_aligned,  fullfile(fig_dir,['group_coherence_spectra_global' cfg.fig_suffix '.png']));
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
    peak_coh_table_full = array2table(...
        [peak_coh_brainEMG, peak_coh_brainSpine, peak_coh_spineEMG, ...
         peak_coh_brainEMG_nt, peak_coh_brainSpine_nt, peak_coh_spineEMG_nt], ...
        'VariableNames',{'BrainEMG_MT','BrainSpine_MT','SpineEMG_MT', ...
        'BrainEMG_NT','BrainSpine_NT','SpineEMG_NT'},'RowNames',subs);
    coh_summary = [median(peak_coh_brainEMG,'omitnan'),    median(peak_coh_brainSpine,'omitnan'),    median(peak_coh_spineEMG,'omitnan'), ...
        median(peak_coh_brainEMG_nt,'omitnan'), median(peak_coh_brainSpine_nt,'omitnan'), median(peak_coh_spineEMG_nt,'omitnan'); ...
        mad(peak_coh_brainEMG,1),    mad(peak_coh_brainSpine,1),    mad(peak_coh_spineEMG,1), ...
        mad(peak_coh_brainEMG_nt,1), mad(peak_coh_brainSpine_nt,1), mad(peak_coh_spineEMG_nt,1)];
    supp_table_1b = [peak_coh_table_full; array2table(coh_summary, ...
        'VariableNames',{'BrainEMG_MT','BrainSpine_MT','SpineEMG_MT', ...
        'BrainEMG_NT','BrainSpine_NT','SpineEMG_NT'},'RowNames',{'Median','MAD'})];
    fprintf('\n  Supplementary Table 1b: Peak coherence (R2, lag mode: %s)\n', lag_mode); disp(supp_table_1b);
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
    
    %% Directed coherence bar plot
    forward_areas = zeros(nSubs,3); reverse_areas = zeros(nSubs,3);
    for ss = 1:nSubs
        forward_areas(ss,:) = [results(ss).brainEMG.forward_area, results(ss).brainSpine.forward_area, results(ss).spineEMG.forward_area];
        reverse_areas(ss,:) = [results(ss).brainEMG.reverse_area, results(ss).brainSpine.reverse_area, results(ss).spineEMG.reverse_area];
    end
    mean_fwd = mean(forward_areas,1); mean_rev = mean(reverse_areas,1);
    sem_fwd  = std(forward_areas,0,1)/sqrt(nSubs); sem_rev = std(reverse_areas,0,1)/sqrt(nSubs);
    data_bar = [mean_fwd; mean_rev]'; sems_bar = [sem_fwd; sem_rev]';
    x_labels = {'Brain\rightarrowEMG','Brain\rightarrowSpine','Spine\rightarrowEMG'};
    hfig_bar = figure('Color','w'); hold on;
    b = bar(data_bar,'grouped');
    b(1).FaceColor = [0.3 0.6 0.9]; b(1).FaceAlpha = 0.5;
    b(2).FaceColor = [0.9 0.3 0.3]; b(2).FaceAlpha = 0.5;
    [ngroups, nbars] = size(data_bar);
    groupwidth = min(0.8, nbars/(nbars+1.5));
    jitter = 0.015;
    for i = 1:nbars
        xpos = (1:ngroups) - groupwidth/2 + (2*i-1)*groupwidth/(2*nbars);
        errorbar(xpos, data_bar(:,i), sems_bar(:,i), 'k.','LineWidth',1.5);
        for j = 1:ngroups
            y = forward_areas(:,j); if i==2, y=reverse_areas(:,j); end
            xj = xpos(j) + jitter*randn(nSubs,1);
            scatter(xj, y, 45,'MarkerFaceColor',[0.35 0.35 0.35],'MarkerEdgeColor','none','MarkerFaceAlpha',0.55);
            scatter(xj(p1_idx), y(p1_idx), 180,'o','MarkerEdgeColor','k','LineWidth',2,'MarkerFaceColor','none');
        end
    end
    set(gca,'XTick',1:ngroups,'XTickLabel',x_labels,'FontSize',14);
    ylabel('Coherence area (10-35 Hz)','FontSize',14);
    ylim([0 max(data_bar(:))*1.3]);
    legend({'Forward','Reverse'},'Location','northwest','FontSize',14); box off;
    title('Group directed coherence (BS)','Interpreter','none','FontSize',14);
    if saveFigs
        savefig(hfig_bar, fullfile(fig_dir,['group_directed_coherence' cfg.fig_suffix '.fig']));
    end
    
    %% Directionality comparison
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
    
    end  % run_coherence_spectra
    
    
    %% =========================================================================
    %  UTILITY FUNCTIONS
    %% =========================================================================
    
    function coh_vec = get_coh(fmat, lag_mode)
    % Extract "the" coherence vector from a sp2a2_R2_mt output matrix
    % according to lag_mode.
    %   'full'    - column 4, the standard total coherence at all lags
    %               including zero.
    %   'nonzero' - columns 11+12 summed, the forward+reverse decomposition
    %               from Halliday (2015), i.e. coherence excluding the
    %               instantaneous/zero-lag component.
    % These two (plus the zero-lag-only column, unused here) sum back to
    % the same total coherence.
    switch lag_mode
        case 'full'
            coh_vec = fmat(:,4);
        case 'nonzero'
            coh_vec = fmat(:,11) + fmat(:,12);
        otherwise
            error('Unknown lag_mode: %s (expected ''full'' or ''nonzero'')', lag_mode);
    end
    end
    
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
    
    function thresh_vec = compute_surrogate_threshold(x_ts, y_ts, samp_rate, seg_pwr, n_surr, opt_str, pctile_val, lag_mode)
    % Pointwise surrogate threshold at each frequency bin, built from the
    % SAME coherence quantity (lag_mode) being tested against it, so the
    % null distribution actually matches the statistic under test.
    % pctile_val: percentile (default 95). Significance = exceeds threshold at
    % ANY frequency in the band (no multiple comparisons correction).
    if nargin < 7 || isempty(pctile_val), pctile_val = 95; end
    if nargin < 8 || isempty(lag_mode), lag_mode = 'nonzero'; end
    surr_coh_mat = [];
    n_samp = length(y_ts);
    Y_fft  = fft(y_ts);
    for si = 1:n_surr
        Y_rand_fft = Y_fft;
        if mod(n_samp, 2) == 0
            n_pos  = n_samp/2 - 1;
            phases = exp(1i * 2*pi * rand(1, n_pos));
            Y_rand_fft(2 : n_samp/2)     = abs(Y_fft(2 : n_samp/2)) .* phases;
            Y_rand_fft(n_samp/2+2 : end) = conj(fliplr(Y_rand_fft(2 : n_samp/2)));
        else
            n_pos  = (n_samp-1)/2;
            phases = exp(1i * 2*pi * rand(1, n_pos));
            Y_rand_fft(2 : (n_samp+1)/2)   = abs(Y_fft(2 : (n_samp+1)/2)) .* phases;
            Y_rand_fft((n_samp+3)/2 : end) = conj(fliplr(Y_rand_fft(2 : (n_samp+1)/2)));
        end
        Y_rand = real(ifft(Y_rand_fft));
        try
            [f_surr, ~, ~] = sp2a2_R2_mt(x_ts', Y_rand', samp_rate, seg_pwr, opt_str);
            coh_surr = get_coh(f_surr, lag_mode)';
            if isempty(surr_coh_mat)
                surr_coh_mat = nan(n_surr, length(coh_surr));
            end
            surr_coh_mat(si, :) = coh_surr;
        catch
            continue
        end
    end
    thresh_vec = prctile(surr_coh_mat, pctile_val, 1);
    end