%% RUN_SPECTRA_BS.m
    % Downstream spectra and coherence analysis — BSLaw forward model.
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
    cfg.subs_coh   = intersect(cfg.subs_brain, cfg.subs_spine);
    
    ve_configs = {
        'bslaw_prevalence', 'VE_spine_prevalence_sub%s_forspectra_BS.mat', 'VE';
        };
    
    cfg.brain_ve_suffix  = '_brain_pct10';
    cfg.spine_ve_suffix  = '_BS';
    cfg.cluster_method   = 'manual';
    cfg.manual_boundary  = 25;
    
    cfg.run_surrogates  = false;    % set false to skip (uses NeuroSpec CL as fallback)
    cfg.n_surrogates    = 500;
    cfg.surr_pctile_BE  = 95;
    cfg.surr_pctile_BS  = 95;     
    cfg.surr_pctile_SE  = 95;
    
    cfg.fband         = [10 35];   % analysis band (significance, power, coherence, PSI)
    cfg.cluster_fband = [10 40];   % wider band for peak extraction / GMM clustering only
    
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
    fprintf('  Cluster method: %s', cfg.cluster_method)
    if strcmp(cfg.cluster_method,'manual')
        fprintf(' (boundary=%.0f Hz)', cfg.manual_boundary)
    end
    fprintf('\n  Analysis band:  %.0f-%.0f Hz', cfg.fband(1), cfg.fband(2))
    fprintf('\n  Cluster band:   %.0f-%.0f Hz', cfg.cluster_fband(1), cfg.cluster_fband(2))
    fprintf('\n  Run surrogates: %d', cfg.run_surrogates)
    fprintf('\n  Surr percentiles: BE=%d  BS=%d  SE=%d\n', ...
        cfg.surr_pctile_BE, cfg.surr_pctile_BS, cfg.surr_pctile_SE)
    fprintf('===========================\n\n')
    
    fprintf('\n>>> STEP B: Coherence spectra\n')
    for vi = 1:size(ve_configs,1)
        cfg.spine_ve_pattern = ve_configs{vi,2};
        cfg.spine_ve_varname = ve_configs{vi,3};
        cfg.fig_suffix       = ['_' ve_configs{vi,1}];
        cfg.fig_dir = fullfile(cfg.save_dir, 'figures', ['spectra_' ve_configs{vi,1}]);
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
    cluster_fband   = cfg.cluster_fband;
    seg_pwr         = cfg.seg_pwr;
    lat_min         = cfg.lat_min;
    lat_max         = cfg.lat_max;
    saveFigs        = cfg.saveFigs;
    fig_dir         = cfg.fig_dir;
    cluster_method  = cfg.cluster_method;
    manual_boundary = cfg.manual_boundary;
    n_surrogates    = cfg.n_surrogates;
    run_surrogates  = cfg.run_surrogates;
    surr_pctile_BE  = cfg.surr_pctile_BE;
    surr_pctile_BS  = cfg.surr_pctile_BS;
    surr_pctile_SE  = cfg.surr_pctile_SE;
    
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
                statBcont, statEMGcont, samp_rate, seg_pwr, n_surrogates, opt_str, surr_pctile_BE);
            results(ss).surr_thresh_BS = compute_surrogate_threshold(...
                statBcont, statScont,   samp_rate, seg_pwr, n_surrogates, opt_str, surr_pctile_BS);
            results(ss).surr_thresh_SE = compute_surrogate_threshold(...
                statScont, statEMGcont, samp_rate, seg_pwr, n_surrogates, opt_str, surr_pctile_SE);
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
            coh_BE  = (f1_mt(:,11)+f1_mt(:,12))*1e3;
            coh_BS  = (f2_mt(:,11)+f2_mt(:,12))*1e3;
            coh_SE  = (f3_mt(:,11)+f3_mt(:,12))*1e3;
            pairs      = {coh_BE, coh_BS, coh_SE};
            cols       = {col_BE, col_BS, col_SE};
            %linestyles = {'-','--',':'};
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
    %             plot(fax_p1, surr_thresh_p1{pp}, '--', 'LineWidth',1, ...
    %                 'Color',[0 0 0 0.5], 'HandleVisibility','off');
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
            all_dir_brainEMG_stat   = nan(nSubs, nFreqs);
            all_dir_brainSpine_stat = nan(nSubs, nFreqs);
            all_dir_spineEMG_stat   = nan(nSubs, nFreqs);
        end
    
        freq_axis = f1_mt(:,1)';
        all_dir_brainEMG_stat(ss,:)   = (f1_mt(:,11) + f1_mt(:,12))';
        all_dir_brainSpine_stat(ss,:) = (f2_mt(:,11) + f2_mt(:,12))';
        all_dir_spineEMG_stat(ss,:)   = (f3_mt(:,11) + f3_mt(:,12))';
    
        %% Peak coherence / frequency — analysis band (for supplementary table)
        band_mask    = freq_axis >= fband(1) & freq_axis <= fband(2);
        freqs_inband = freq_axis(band_mask);
    
        [pk1, i1] = max(f1_mt(band_mask, 4));
        peak_freq_brainEMG(ss)   = freqs_inband(i1); peak_coh_brainEMG(ss)   = pk1;
        [pk2, i2] = max(f2_mt(band_mask, 4));
        peak_freq_brainSpine(ss) = freqs_inband(i2); peak_coh_brainSpine(ss) = pk2;
        [pk3, i3] = max(f3_mt(band_mask, 4));
        peak_freq_spineEMG(ss)   = freqs_inband(i3); peak_coh_spineEMG(ss)   = pk3;
    
        %% Peak frequencies — cluster band (wider, for GMM only)
        clust_mask  = freq_axis >= cluster_fband(1) & freq_axis <= cluster_fband(2);
        freqs_clust = freq_axis(clust_mask);
        [~, ic1] = max(f1_mt(clust_mask, 4));
        [~, ic2] = max(f2_mt(clust_mask, 4));
        [~, ic3] = max(f3_mt(clust_mask, 4));
        results(ss).clust_peak_BE = freqs_clust(ic1);
        results(ss).clust_peak_BS = freqs_clust(ic2);
        results(ss).clust_peak_SE = freqs_clust(ic3);
    
        results(ss).dir_BE = (f1_mt(:,11) + f1_mt(:,12))';
        results(ss).dir_BS = (f2_mt(:,11) + f2_mt(:,12))';
        results(ss).dir_SE = (f3_mt(:,11) + f3_mt(:,12))';
    
        %% NT peak coherence (analysis band)
        fax_nt_pk       = f1_nt(:,1)';
        band_nt         = fax_nt_pk >= fband(1) & fax_nt_pk <= fband(2);
        freqs_nt_inband = fax_nt_pk(band_nt);
        peak_coh_brainEMG_nt(ss)   = max(f1_nt(band_nt, 4));
        peak_coh_brainSpine_nt(ss) = max(f2_nt(band_nt, 4));
        peak_coh_spineEMG_nt(ss)   = max(f3_nt(band_nt, 4));
        [~, i1_nt] = max(f1_nt(band_nt, 4));
        [~, i2_nt] = max(f2_nt(band_nt, 4));
        [~, i3_nt] = max(f3_nt(band_nt, 4));
        results(ss).peak_freq_BE_nt = freqs_nt_inband(i1_nt);
        results(ss).peak_freq_BS_nt = freqs_nt_inband(i2_nt);
        results(ss).peak_freq_SE_nt = freqs_nt_inband(i3_nt);
    
        results(ss).sub                 = sub;
        results(ss).partial_mt_BE_spine = fp_mt_BE_spine';
        results(ss).full_coh_mt_BE      = fc_mt_BE';
        results(ss).full_coh_mt_SE      = fc_mt_SE';
        results(ss).full_coh_mt_BS      = fc_mt_BS';
        results(ss).freq_axis_mt        = f1_mt(:,1)';
        results(ss).peak_freq_BE_mt     = peak_freq_brainEMG(ss);
        results(ss).conf95_BE           = cl1_mt.ch_c95;
        results(ss).conf95_BS           = cl2_mt.ch_c95;
        results(ss).conf95_SE           = cl3_mt.ch_c95;
    
        %% Directed coherence areas
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
    %  CLUSTERING  (uses cluster_fband peaks, not fband peaks)
    %% =========================================================================
    p1_idx = 1;
    
    clust_peak_BE_all = arrayfun(@(s) s.clust_peak_BE, results)';
    clust_peak_BS_all = arrayfun(@(s) s.clust_peak_BS, results)';
    clust_peak_SE_all = arrayfun(@(s) s.clust_peak_SE, results)';
    all_peaks_pooled  = [clust_peak_BE_all; clust_peak_BS_all; clust_peak_SE_all];
    all_peaks_pooled  = all_peaks_pooled(~isnan(all_peaks_pooled));
    
    fprintf('\n=== GMM clustering of peak frequencies (cluster band: %.0f-%.0f Hz) ===\n', ...
        cluster_fband(1), cluster_fband(2));
    fprintf('  Pooled peaks (n=%d): ', numel(all_peaks_pooled));
    fprintf('%.1f  ', sort(all_peaks_pooled)'); fprintf('\n');
    
    max_k = 4; bic = nan(1,max_k); gm_fit = cell(1,max_k);
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
    gm_best = gm_fit{best_k};
    [~, sort_idx] = sort(gm_best.mu);
    cluster_means = gm_best.mu(sort_idx);
    fprintf('  BIC: '); fprintf('k=%d: %.1f  ',[1:max_k; bic]); fprintf('\n');
    fprintf('  Best k = %d\n', best_k);
    
    boundaries   = [cluster_fband(1), arrayfun(@(i) mean(cluster_means(i:i+1)), 1:best_k-1), cluster_fband(2)];
    nClusters    = best_k;
    cluster_cols = lines(nClusters);
    
    if nClusters < 2
        warning('GMM k=1; falling back to largest-gap split.');
        sp_ = sort(all_peaks_pooled); [~,gi] = max(diff(sp_));
        boundaries    = [cluster_fband(1), mean(sp_(gi:gi+1)), cluster_fband(2)];
        nClusters     = 2;
        cluster_means = [mean(sp_(sp_<boundaries(2))), mean(sp_(sp_>=boundaries(2)))];
        cluster_cols  = lines(2);
    end
    
    if strcmp(cluster_method, 'manual')
        nClusters     = 2;
        boundaries    = [cluster_fband(1), manual_boundary, cluster_fband(2)];
        cluster_means = [mean([cluster_fband(1) manual_boundary]), mean([manual_boundary cluster_fband(2)])];
        cluster_cols  = lines(2);
        fprintf('  Using manual boundaries: [%.0f %.0f %.0f] Hz\n', boundaries);
    end
    fprintf('  Boundaries (Hz): '); fprintf('%.1f  ', boundaries); fprintf('\n');
    fprintf('  Cluster means (Hz): '); fprintf('%.1f  ', cluster_means); fprintf('\n');
    fprintf('======================================\n\n');
    
    %% Sub-band peak frequencies per pair (using cluster_fband peaks per results)
    peak_freq_clust = nan(nSubs, 3, nClusters);
    for ss = 1:nSubs
        fax = results(ss).freq_axis_mt;
        if isempty(fax), continue; end
        pairs_dir = {results(ss).dir_BE, results(ss).dir_BS, results(ss).dir_SE};
        for p = 1:3
            d = pairs_dir{p};
            if isempty(d) || all(isnan(d)), continue; end
            for c = 1:nClusters
                mask_c = fax >= boundaries(c) & fax <= boundaries(c+1);
                f_c = fax(mask_c);
                if sum(mask_c) < 2, continue; end
                [~, ic] = max(d(mask_c));
                peak_freq_clust(ss,p,c) = f_c(ic);
            end
        end
    end
    
    pair_names_r = {'Brain-EMG','Brain-Spine','Spine-EMG'};
    fprintf('  Per-pair peak frequencies by cluster:\n');
    for c = 1:nClusters
        fprintf('  Cluster %d (%.1f-%.1f Hz, mean=%.1f Hz):\n', ...
            c, boundaries(c), boundaries(c+1), cluster_means(c));
        for p = 1:3
            pf = peak_freq_clust(:,p,c); pf = pf(~isnan(pf));
            if isempty(pf), fprintf('    %s: no data\n', pair_names_r{p}); continue; end
            fprintf('    %s: median=%.1f  range=[%.1f %.1f]  n=%d\n', ...
                pair_names_r{p}, median(pf), min(pf), max(pf), numel(pf));
        end
    end
    
    %% Within-subject frequency differences
    fprintf('\n=== Within-subject peak frequency differences ===\n');
    for c = 1:nClusters
        fprintf('\n  Cluster %d (%.0f-%.0f Hz):\n', c, boundaries(c), boundaries(c+1));
        diff_BE_BS = abs(peak_freq_clust(:,1,c) - peak_freq_clust(:,2,c));
        diff_BE_SE = abs(peak_freq_clust(:,1,c) - peak_freq_clust(:,3,c));
        diff_BS_SE = abs(peak_freq_clust(:,2,c) - peak_freq_clust(:,3,c));
        ok_BS = ~isnan(diff_BE_BS); ok_SE = ~isnan(diff_BE_SE); ok_SS = ~isnan(diff_BS_SE);
        half_bw = (boundaries(c+1) - boundaries(c)) / 2;
        if sum(ok_BS)>0
            fprintf('  |BE-BS|: median=%.2f  IQR=[%.2f %.2f]  n=%d  pct<half_bw=%.0f%%\n', ...
                median(diff_BE_BS(ok_BS)), prctile(diff_BE_BS(ok_BS),25), ...
                prctile(diff_BE_BS(ok_BS),75), sum(ok_BS), 100*mean(diff_BE_BS(ok_BS)<half_bw));
        end
        if sum(ok_SE)>0
            fprintf('  |BE-SE|: median=%.2f  IQR=[%.2f %.2f]  n=%d  pct<half_bw=%.0f%%\n', ...
                median(diff_BE_SE(ok_SE)), prctile(diff_BE_SE(ok_SE),25), ...
                prctile(diff_BE_SE(ok_SE),75), sum(ok_SE), 100*mean(diff_BE_SE(ok_SE)<half_bw));
        end
        if sum(ok_SS)>0
            fprintf('  |BS-SE|: median=%.2f  IQR=[%.2f %.2f]  n=%d  pct<half_bw=%.0f%%\n', ...
                median(diff_BS_SE(ok_SS)), prctile(diff_BS_SE(ok_SS),25), ...
                prctile(diff_BS_SE(ok_SS),75), sum(ok_SS), 100*mean(diff_BS_SE(ok_SS)<half_bw));
        end
    end
    
    %% Subject cluster classification
    subj_cluster = nan(nSubs, 1);
    for ss = 1:nSubs
        fax = results(ss).freq_axis_mt;
        if isempty(fax), continue; end
        d = results(ss).dir_BE;
        if isempty(d) || all(isnan(d)), continue; end
        best_coh = -inf;
        for c = 1:nClusters
            mask_c = fax >= boundaries(c) & fax <= boundaries(c+1);
            pk = max(d(mask_c));
            if pk > best_coh
                best_coh = pk;
                subj_cluster(ss) = c;
            end
        end
    end
    
    %% Cluster consistency across pairs
    fprintf('\n=== Cluster consistency across signal pairs ===\n');
    fprintf('  Sub           BE    BS    SE    Consistent\n');
    n_consistent = 0;
    for ss = 1:nSubs
        fax = results(ss).freq_axis_mt;
        pairs_dir = {results(ss).dir_BE, results(ss).dir_BS, results(ss).dir_SE};
        pair_clusters = nan(1,3);
        for p = 1:3
            d = pairs_dir{p};
            if isempty(d) || all(isnan(d)), continue; end
            best_coh = -inf;
            for c = 1:nClusters
                mask_c = fax >= boundaries(c) & fax <= boundaries(c+1);
                pk = max(d(mask_c));
                if pk > best_coh
                    best_coh = pk;
                    pair_clusters(p) = c;
                end
            end
        end
        c_BE = pair_clusters(1); c_BS = pair_clusters(2); c_SE = pair_clusters(3);
        consistent = ~isnan(c_BE) && ~isnan(c_BS) && ~isnan(c_SE) && c_BE==c_BS && c_BE==c_SE;
        if consistent, n_consistent = n_consistent + 1; end
        fprintf('  %-12s  C%d    C%d    C%d    %s\n', subs{ss}, ...
            c_BE, c_BS, c_SE, conditional_str(consistent,'YES','no'));
    end
    fprintf('  Consistent: %d/%d subjects\n', n_consistent, nSubs);
    fprintf('================================================\n');
    
    %% =========================================================================
    %  FIGURES
    %% =========================================================================
    prominence_factor = 1.5;
    tolerance_hz      = 2;
    
    coloc_BS = nan(nSubs, nClusters);
    coloc_SE = nan(nSubs, nClusters);
    dist_BS  = nan(nSubs, nClusters);
    dist_SE  = nan(nSubs, nClusters);
    sig_BE   = nan(nSubs, 1);
    sig_BS   = nan(nSubs, 1);
    sig_SE   = nan(nSubs, 1);
    
    for ss = 1:nSubs
        fax = results(ss).freq_axis_mt;
        if isempty(fax), continue; end
        band_mask = fax >= fband(1) & fax <= fband(2);
        sig_BE(ss) = any(results(ss).dir_BE(band_mask) > results(ss).surr_thresh_BE(band_mask));
        sig_BS(ss) = any(results(ss).dir_BS(band_mask) > results(ss).surr_thresh_BS(band_mask));
        sig_SE(ss) = any(results(ss).dir_SE(band_mask) > results(ss).surr_thresh_SE(band_mask));
    end
    
    fprintf('\n=== Beta band coherence significance (surrogate threshold) ===\n');
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
    
    %% KDE + scatter  (axes use cluster_fband)
    kde_bw = 1.5;
    f_eval = linspace(cluster_fband(1)-2, cluster_fband(2)+2, 300);
    hfig_clust = figure('Color','w','Position',[100 100 1100 600]);
    markers = {'o','d','s','^'};
    scatter_cfg = {[1 2],'Brain-EMG vs Brain-Spine','Brain-EMG peak (Hz)','Brain-Spine peak (Hz)'; ...
        [1 3],'Brain-EMG vs Spine-EMG',  'Brain-EMG peak (Hz)','Spine-EMG peak (Hz)'; ...
        [2 3],'Brain-Spine vs Spine-EMG','Brain-Spine peak (Hz)','Spine-EMG peak (Hz)'};
    
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
                plot([pf_c(ss) pf_c(ss)],[0 rug_y],'-','Color',cluster_cols(c,:),'LineWidth',lw_r);
            end
        end
        for bi = 2:nClusters
            xline(boundaries(bi),'--k','Alpha',0.3,'HandleVisibility','off');
        end
        xlim([cluster_fband(1)-1 cluster_fband(2)+1]);
        xlabel('Peak frequency (Hz)','FontSize',11); ylabel('Density','FontSize',11);
        title(pair_names_r{p},'FontSize',12); set(gca,'FontSize',11); box off;
        if p==1
            leg_e = arrayfun(@(c) sprintf('Cluster %d (%.0f-%.0f Hz)', ...
                c,boundaries(c),boundaries(c+1)), 1:nClusters,'UniformOutput',false);
            legend(leg_e,'Location','northwest','Box','off','FontSize',9);
        end
    end
    
    dlim = [cluster_fband(1)-1 cluster_fband(2)+1];
    for sp = 1:3
        subplot(2,4,4+sp); hold on;
        pi1 = scatter_cfg{sp,1}(1); pi2 = scatter_cfg{sp,1}(2);
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
                plot(xv(ok_xy),yv(ok_xy),'-','Color',[0.8 0.8 0.8],'LineWidth',0.8,'HandleVisibility','off');
            end
            for c = 1:nClusters
                xc = peak_freq_clust(ss,pi1,c); yc = peak_freq_clust(ss,pi2,c);
                if isnan(xc)||isnan(yc), continue; end
                scatter(xc,yc,ms,cluster_cols(c,:),markers{c},'filled', ...
                    'MarkerEdgeColor','w','LineWidth',lw,'MarkerFaceAlpha',0.75,'HandleVisibility','off');
            end
        end
        hleg = gobjects(nClusters,1);
        for c = 1:nClusters
            hleg(c) = scatter(nan,nan,60,cluster_cols(c,:),markers{c},'filled','MarkerFaceAlpha',0.75);
        end
        legend(hleg, arrayfun(@(c) sprintf('Cluster %d',c),1:nClusters,'UniformOutput',false), ...
            'Location','northwest','Box','off','FontSize',10);
        xlim(dlim); ylim(dlim); axis square;
        xlabel(scatter_cfg{sp,3},'FontSize',11); ylabel(scatter_cfg{sp,4},'FontSize',11);
        title(scatter_cfg{sp,2},'FontSize',12); set(gca,'FontSize',11); box off;
    end
    
    subplot(2,4,8); axis off;
    ann = [{'Clusters:'}, ...
        arrayfun(@(c) sprintf('C%d: %.0f-%.0f Hz (mean %.1f)', ...
        c,boundaries(c),boundaries(c+1),cluster_means(c)), ...
        1:nClusters,'UniformOutput',false), ...
        {sprintf('KDE bw=%.1f Hz | Large=P1', kde_bw)}, ...
        {sprintf('Cluster band: %.0f-%.0f Hz', cluster_fband(1), cluster_fband(2))}];
    text(0.5,0.55,ann,'HorizontalAlignment','center','FontSize',10,'Units','normalized');
    sgtitle('Peak frequency clustering (BS)','FontSize',13);
    if saveFigs
        savefig(hfig_clust, fullfile(fig_dir,['group_peak_frequency_clustering' cfg.fig_suffix '.fig']));
        saveas(hfig_clust,  fullfile(fig_dir,['group_peak_frequency_clustering' cfg.fig_suffix '.png']));
    end
    
    %% Individual subject spectra — one figure per cluster
    for c = 1:nClusters
        hfig_indiv = figure('Color','w','Position',[100 100 nSubs*180 420]);
        for ss = 1:nSubs
            subplot(1,nSubs,ss); hold on;
            fax = results(ss).freq_axis_mt;
            if isempty(fax), continue; end
            pairs_data  = {results(ss).dir_BE, results(ss).dir_BS, results(ss).dir_SE};
            pair_cols_i = {col_BE, col_BS, col_SE};
            pair_ls_i   = {'-','--',':'};
            band_lo     = boundaries(c);
            band_hi     = boundaries(c+1);
            band_mask_c = fax >= band_lo & fax <= band_hi;
            f_band_c    = fax(band_mask_c);
    
            d_BE         = pairs_data{1};
            BE_peak_freq = nan;
            if ~isempty(d_BE) && ~all(isnan(d_BE)) && sum(band_mask_c) >= 2
                d_BE_band = d_BE(band_mask_c);
                [pks_be, locs_be] = findpeaks(d_BE_band, f_band_c, ...
                    'MinPeakProminence', prominence_factor * median(d_BE_band));
                if ~isempty(pks_be)
                    [~, imax] = max(pks_be); BE_peak_freq = locs_be(imax);
                else
                    [~, imax] = max(d_BE_band); BE_peak_freq = f_band_c(imax);
                end
            end
    
            for p = 1:3
                d = pairs_data{p};
                if isempty(d) || all(isnan(d)), continue; end
                plot(fax, d*1e3, pair_ls_i{p}, ...
                    'Color',[pair_cols_i{p} 0.25],'LineWidth',0.8,'HandleVisibility','off');
                if sum(band_mask_c) >= 2
                    plot(f_band_c, d(band_mask_c)*1e3, pair_ls_i{p}, ...
                        'Color',pair_cols_i{p},'LineWidth',2,'HandleVisibility','off');
                end
                if sum(band_mask_c) < 2, continue; end
                d_band = d(band_mask_c);
                if p == 1
                    [pks, locs_f] = findpeaks(d_band, f_band_c, ...
                        'MinPeakProminence', prominence_factor * median(d_band));
                    if ~isempty(pks)
                        [~,imax] = max(pks);
                        scatter(locs_f(imax), pks(imax)*1e3, 100, pair_cols_i{p}, ...
                            'filled','MarkerEdgeColor','k','LineWidth',1,'HandleVisibility','off');
                        plot([locs_f(imax) locs_f(imax)],[0 pks(imax)*1e3*1.1],':', ...
                            'Color',pair_cols_i{p},'LineWidth',0.8,'HandleVisibility','off');
                    else
                        [pk_val,pk_i] = max(d_band);
                        scatter(f_band_c(pk_i), pk_val*1e3, 80, pair_cols_i{p}, ...
                            'MarkerEdgeColor',pair_cols_i{p},'LineWidth',1.5,'HandleVisibility','off');
                    end
                else
                    [pks, locs_f] = findpeaks(d_band, f_band_c, ...
                        'MinPeakProminence', prominence_factor * median(d_band));
                    if ~isnan(BE_peak_freq) && ~isempty(pks)
                        in_window = abs(locs_f - BE_peak_freq) <= tolerance_hz;
                        if any(in_window)
                            [~,ibest] = max(pks(in_window));
                            locs_win = locs_f(in_window); pks_win = pks(in_window);
                            best_freq = locs_win(ibest); best_val = pks_win(ibest);
                            scatter(best_freq, best_val*1e3, 100, pair_cols_i{p}, ...
                                'filled','MarkerEdgeColor','k','LineWidth',1,'HandleVisibility','off');
                            plot([best_freq best_freq],[0 best_val*1e3*1.1],':', ...
                                'Color',pair_cols_i{p},'LineWidth',0.8,'HandleVisibility','off');
                            if p==2, coloc_BS(ss,c)=1; dist_BS(ss,c)=abs(best_freq-BE_peak_freq);
                            else,    coloc_SE(ss,c)=1; dist_SE(ss,c)=abs(best_freq-BE_peak_freq); end
                        else
                            [~,imax] = max(pks);
                            scatter(locs_f(imax), pks(imax)*1e3, 80, pair_cols_i{p}, ...
                                'MarkerEdgeColor',pair_cols_i{p},'LineWidth',1.5,'HandleVisibility','off');
                            if p==2, coloc_BS(ss,c)=0; dist_BS(ss,c)=min(abs(locs_f-BE_peak_freq));
                            else,    coloc_SE(ss,c)=0; dist_SE(ss,c)=min(abs(locs_f-BE_peak_freq)); end
                        end
                    elseif ~isnan(BE_peak_freq)
                        [pk_val,pk_i] = max(d_band);
                        scatter(f_band_c(pk_i), pk_val*1e3, 80, pair_cols_i{p}, ...
                            'MarkerEdgeColor',pair_cols_i{p},'LineWidth',1.5,'HandleVisibility','off');
                        if p==2, coloc_BS(ss,c)=0; else, coloc_SE(ss,c)=0; end
                    end
                end
            end
    
            yl = ylim; if yl(2)<=yl(1), yl(2)=yl(1)+0.1; end
            fill([band_lo band_hi band_hi band_lo],[yl(1) yl(1) yl(2) yl(2)], ...
                cluster_cols(c,:),'FaceAlpha',0.08,'EdgeColor','none','HandleVisibility','off');
            xline(band_lo,'--','Color',[0.5 0.5 0.5],'Alpha',0.4,'HandleVisibility','off');
            xline(band_hi,'--','Color',[0.5 0.5 0.5],'Alpha',0.4,'HandleVisibility','off');
            xlim([cluster_fband(1)-2 cluster_fband(2)+2]);
            xlabel('Hz','FontSize',9);
            if ss==1
                ylabel('Coherence (x1e-3)','FontSize',10);
                plot(nan,nan,'-','Color',col_BE,'LineWidth',2,'DisplayName','Brain-EMG');
                plot(nan,nan,'--','Color',col_BS,'LineWidth',2,'DisplayName','Brain-Spine');
                plot(nan,nan,':','Color',col_SE,'LineWidth',2,'DisplayName','Spine-EMG');
                plot(nan,nan,'o','MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor','k', ...
                    'MarkerSize',8,'DisplayName','Co-located peak');
                plot(nan,nan,'o','MarkerFaceColor','none','MarkerEdgeColor',[0.5 0.5 0.5], ...
                    'MarkerSize',8,'DisplayName','No co-located peak');
                legend('Location','northeast','Box','off','FontSize',7);
            end
            bs_str = ''; se_str = '';
            if ~isnan(coloc_BS(ss,c)), bs_str = conditional_str(coloc_BS(ss,c)==1,'+','-'); end
            if ~isnan(coloc_SE(ss,c)), se_str = conditional_str(coloc_SE(ss,c)==1,'+','-'); end
            bs_sig_str = conditional_str(sig_BS(ss)==1,'*',' ');
            se_sig_str = conditional_str(sig_SE(ss)==1,'*',' ');
            if ss==p1_idx
                title(sprintf('P%d*\nBS:%s%s SE:%s%s',ss,bs_sig_str,bs_str,se_sig_str,se_str), ...
                    'FontSize',9,'FontWeight','bold');
            else
                title(sprintf('P%d\nBS:%s%s SE:%s%s',ss,bs_sig_str,bs_str,se_sig_str,se_str),'FontSize',9);
            end
            set(gca,'FontSize',9); box off;
        end
        sgtitle(sprintf('Individual spectra — Cluster %d (%.0f-%.0f Hz)\nFilled=co-located (+-%.0f Hz), Hollow=no co-located peak', ...
            c,boundaries(c),boundaries(c+1),tolerance_hz),'FontSize',11,'Color',cluster_cols(c,:));
        if saveFigs
            savefig(hfig_indiv, fullfile(fig_dir,sprintf('individual_spectra_cluster%d%s.fig',c,cfg.fig_suffix)));
            saveas(hfig_indiv,  fullfile(fig_dir,sprintf('individual_spectra_cluster%d%s.png',c,cfg.fig_suffix)));
        end
    end
    
    %% Co-location summary
    fprintf('\n=== Peak co-location summary (+-%.0f Hz of Brain-EMG peak) ===\n', tolerance_hz);
    for c = 1:nClusters
        fprintf('\n  Cluster %d (%.0f-%.0f Hz):\n', c, boundaries(c), boundaries(c+1));
        fprintf('  Sub           BS_coloc  BS_dist(Hz)  SE_coloc  SE_dist(Hz)\n');
        for ss = 1:nSubs
            fprintf('  %-12s  %s         %s            %s         %s\n', subs{ss}, ...
                fmt_coloc(coloc_BS(ss,c)), fmt_dist(dist_BS(ss,c)), ...
                fmt_coloc(coloc_SE(ss,c)), fmt_dist(dist_SE(ss,c)));
        end
        fprintf('  Brain-Spine co-located: %d/%d (%.0f%%)\n', ...
            sum(coloc_BS(:,c)==1,'omitnan'), sum(~isnan(coloc_BS(:,c))), ...
            100*mean(coloc_BS(:,c)==1,'omitnan'));
        fprintf('  Spine-EMG co-located:   %d/%d (%.0f%%)\n', ...
            sum(coloc_SE(:,c)==1,'omitnan'), sum(~isnan(coloc_SE(:,c))), ...
            100*mean(coloc_SE(:,c)==1,'omitnan'));
        fprintf('  Median BS distance: %.2f Hz\n', median(dist_BS(:,c),'omitnan'));
        fprintf('  Median SE distance: %.2f Hz\n', median(dist_SE(:,c),'omitnan'));
    end
    fprintf('======================================================\n');
    
    %% Grouped spectra figure
    hfig_grouped = figure('Color','w','Position',[100 100 900 nClusters*280]);
    pair_cols_grp  = {col_BE, col_BS, col_SE};
    pair_names_grp = {'Brain-EMG','Brain-Spine','Spine-EMG'};
    for c = 1:nClusters
        subs_in_cluster = find(subj_cluster == c);
        n_in_cluster    = numel(subs_in_cluster);
        for p = 1:3
            ax_idx = (c-1)*3 + p;
            subplot(nClusters,3,ax_idx); hold on;
            col_p = pair_cols_grp{p};
            for si = 1:n_in_cluster
                ss  = subs_in_cluster(si);
                fax = results(ss).freq_axis_mt;
                if isempty(fax), continue; end
                switch p
                    case 1, d = results(ss).dir_BE;
                    case 2, d = results(ss).dir_BS;
                    case 3, d = results(ss).dir_SE;
                end
                if isempty(d) || all(isnan(d)), continue; end
                mask_c = fax >= boundaries(c) & fax <= boundaries(c+1);
                f_c = fax(mask_c); [~,ipk] = max(d(mask_c)); peak_f = f_c(ipk);
                if ss==p1_idx, lw=2.5; lc='k'; else, lw=1.2; lc=col_p; end
                plot(fax, d*1e3, '-','Color',[lc 0.6],'LineWidth',lw);
                plot([peak_f peak_f],[0 max(d(mask_c))*1e3*1.05],':','Color',lc,'LineWidth',0.8,'HandleVisibility','off');
            end
            yl = ylim;
            fill([boundaries(c) boundaries(c+1) boundaries(c+1) boundaries(c)], ...
                [yl(1) yl(1) yl(2) yl(2)],cluster_cols(c,:),'FaceAlpha',0.08,'EdgeColor','none','HandleVisibility','off');
            xline(fband(1),'--k','Alpha',0.2,'HandleVisibility','off');
            xline(fband(2),'--k','Alpha',0.2,'HandleVisibility','off');
            xlim([2 45]);
            xlabel('Frequency (Hz)','FontSize',11); ylabel('Coherence (x1e-3)','FontSize',11);
            if p==2
                title(sprintf('Cluster %d (%.0f-%.0f Hz, n=%d)',c,boundaries(c),boundaries(c+1),n_in_cluster),'FontSize',12);
            else
                title(pair_names_grp{p},'FontSize',12);
            end
            set(gca,'FontSize',11); box off;
        end
    end
    sgtitle('Coherence spectra grouped by beta cluster','FontSize',13);
    if saveFigs
        savefig(hfig_grouped, fullfile(fig_dir,['grouped_spectra_by_cluster' cfg.fig_suffix '.fig']));
        saveas(hfig_grouped,  fullfile(fig_dir,['grouped_spectra_by_cluster' cfg.fig_suffix '.png']));
    end
    
    %% Partial vs full coherence — paired test across participants
    % Same comparison the aligned figure below shows visually (full vs
    % partial Brain-EMG coherence, controlling for spine), but reduced to
    % one paired value per subject: each subject's own full-coherence peak
    % (within fband) vs the partial coherence at that same frequency bin.
    % Both come from mt_partial_coherence (same estimator/tapering), so the
    % comparison isn't confounded by using two different coherence methods.
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
    
    %% Group coherence spectra aligned to peak
    pair_data_g  = {all_dir_brainEMG_stat, all_dir_brainSpine_stat, all_dir_spineEMG_stat};
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
        pairs_dir = {results(ss).dir_BE, results(ss).dir_BS, results(ss).dir_SE};
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
    sgtitle('Group coherence spectra — aligned to each pair''s own beta peak','FontSize',13);
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
    
    %% Supplementary tables
    peak_freq_table_full = array2table(...
        [peak_freq_brainEMG, peak_freq_brainSpine, peak_freq_spineEMG], ...
        'VariableNames',{'BrainEMG_Hz','BrainSpine_Hz','SpineEMG_Hz'},'RowNames',subs);
    freq_summary = [median(peak_freq_brainEMG,'omitnan'), median(peak_freq_brainSpine,'omitnan'), median(peak_freq_spineEMG,'omitnan'); ...
        mad(peak_freq_brainEMG,1), mad(peak_freq_brainSpine,1), mad(peak_freq_spineEMG,1)];
    supp_table_1a = [peak_freq_table_full; array2table(freq_summary, ...
        'VariableNames',{'BrainEMG_Hz','BrainSpine_Hz','SpineEMG_Hz'},'RowNames',{'Median','MAD'})];
    fprintf('\n  Supplementary Table 1a: Peak coherence frequency (Hz, analysis band)\n'); disp(supp_table_1a);
    try
        writetable(supp_table_1a, fullfile(fig_dir,'SuppTable1a_peak_frequency_BS.csv'),'WriteRowNames',true);
    catch ME, warning('Could not save Table 1a: %s', ME.message); end
    
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
    fprintf('\n  Supplementary Table 1b: Peak coherence (R2)\n'); disp(supp_table_1b);
    try
        writetable(supp_table_1b, fullfile(fig_dir,'SuppTable1b_peak_coherence_BS.csv'),'WriteRowNames',true);
    catch ME, warning('Could not save Table 1b: %s', ME.message); end
    
    for c = 1:nClusters
        pf_table = array2table(peak_freq_clust(:,:,c), ...
            'VariableNames',{'BrainEMG_Hz','BrainSpine_Hz','SpineEMG_Hz'},'RowNames',subs);
        pf_summary = [median(peak_freq_clust(:,1,c),'omitnan'), median(peak_freq_clust(:,2,c),'omitnan'), median(peak_freq_clust(:,3,c),'omitnan'); ...
            mad(peak_freq_clust(:,1,c),1), mad(peak_freq_clust(:,2,c),1), mad(peak_freq_clust(:,3,c),1)];
        supp_table_c = [pf_table; array2table(pf_summary, ...
            'VariableNames',{'BrainEMG_Hz','BrainSpine_Hz','SpineEMG_Hz'},'RowNames',{'Median','MAD'})];
        fprintf('\n  Supp Table 1c%d: Peak freq within Cluster %d (%.0f-%.0f Hz)\n', c,c,boundaries(c),boundaries(c+1));
        disp(supp_table_c);
        try
            writetable(supp_table_c, fullfile(fig_dir, ...
                sprintf('SuppTable1c_cluster%d_peak_freq_BS.csv',c)),'WriteRowNames',true);
        catch ME, warning('Could not save Table 1c%d: %s', c, ME.message); end
    end
    
    %% Peak frequency statistics
    freq_mat = [peak_freq_brainEMG, peak_freq_brainSpine, peak_freq_spineEMG];
    complete = all(isfinite(freq_mat), 2);
    freq_mat_complete = freq_mat(complete, :);
    n_complete = sum(complete);
    fprintf('\n=== Peak frequency statistics ===\n');
    fprintf('  N subjects with complete data: %d\n', n_complete);
    pair_names = {'Brain-EMG','Brain-Spine','Spine-EMG'};
    if n_complete >= 3
        [p_fried, tbl_fried, ~] = friedman(freq_mat_complete, 1, 'off');
        fprintf('\n  Friedman test: chi2(%d)=%.3f, p=%.4g\n', tbl_fried{2,3}, tbl_fried{2,2}, p_fried);
        combos = nchoosek(1:3, 2);
        fprintf('\n  Post-hoc pairwise Wilcoxon (uncorrected):\n');
        p_posthoc = nan(size(combos,1),1);
        for c = 1:size(combos,1)
            a = combos(c,1); b = combos(c,2);
            [p_wil, ~] = signrank(freq_mat_complete(:,a), freq_mat_complete(:,b));
            p_posthoc(c) = p_wil;
            fprintf('  %-20s vs %-20s   p=%.4g\n', pair_names{a}, pair_names{b}, p_wil);
        end
        fprintf('  Bonferroni alpha: %.4f\n', 0.05/3);
        for c = 1:size(combos,1)
            a = combos(c,1); b = combos(c,2);
            fprintf('    %s vs %s: p=%.4g [%s]\n', pair_names{a}, pair_names{b}, ...
                p_posthoc(c), conditional_str(p_posthoc(c)<0.05/3,'significant','ns'));
        end
    end
    fprintf('\n  Descriptive (Hz):\n');
    for pp = 1:3
        v = freq_mat(:,pp);
        fprintf('  %-15s  median=%6.2f  MAD=%5.2f  min=%5.2f  max=%5.2f\n', ...
            pair_names{pp}, median(v,'omitnan'), mad(v,1), min(v,[],'omitnan'), max(v,[],'omitnan'));
    end
    fprintf('================================================\n');
    
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
    
    function s = fmt_coloc(v)
    if isnan(v),  s = 'n/a  ';
    elseif v==1,  s = 'YES  ';
    else,         s = 'no   '; end
    end
    
    function s = fmt_dist(v)
    if isnan(v),  s = 'n/a ';
    else,         s = sprintf('%.2f', v); end
    end
    
    function thresh_vec = compute_surrogate_threshold(x_ts, y_ts, samp_rate, seg_pwr, n_surr, opt_str, pctile_val)
    % Pointwise surrogate threshold at each frequency bin.
    % pctile_val: percentile (default 95). Significance = exceeds threshold at
    % ANY frequency in the band (no multiple comparisons correction).
    if nargin < 7 || isempty(pctile_val), pctile_val = 95; end
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
            coh_surr = f_surr(:,4)';
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