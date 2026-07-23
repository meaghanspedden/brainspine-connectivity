%% RUN_UNRECTIFIED_EMG.m
% Sensitivity analysis: repeats DICS coherence pipeline for Participant 1
% (OP00212) using UNRECTIFIED EMG, addressing reviewer comment on the
% rectification debate (McClelland, J Neurosci 2012).
%
% Runs Steps 1 and 2:
%   1. Brain-EMG coherence (singleshell BEM, all channels for leadfield)
%   2. Spine-EMG coherence (BSLaw forward model)
%
% Results saved alongside main pipeline figures with '_unrectified' suffix
% for direct visual comparison with rectified results.
%
% KEY DIFFERENCES FROM MAIN PIPELINE:
%   - EMG is NOT rectified
%   - Brain leadfield uses grad_mm (all channels) not braingrad
%   - Spine uses BSLaw leadfield (leadfield_bs variable)
%   - Channel selection applied before spine DICS (Lf.label + EXG1 only)
%
% =========================================================================

clear all
close all
clc

%% =========================================================================
%  USER CONFIG
%% =========================================================================

% Machine-specific paths live in source/brainspine_config.m — edit that
% file to match your local installation.
repo_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(repo_root, 'source'));
paths = brainspine_config();

cfg.fieldtrip_path  = paths.fieldtrip_path;
cfg.spm_path        = paths.spm_path;
cfg.bsc_source_path = fullfile(repo_root, 'source');

cfg.data_root      = paths.data_root;
cfg.lf_bs_path     = paths.lf_path;
cfg.geomfile       = paths.geomfile;
cfg.T_mat_path     = paths.T_mat_path;
cfg.save_dir       = paths.save_dir;
cfg.fig_dir        = fullfile(cfg.save_dir, 'figures', 'unrectified_emg');

cfg.sub = 'OP00212';

% Analysis parameters — match main pipeline exactly
cfg.fband          = [10 35];
cfg.numpermutation = 500;
cfg.rng_seed       = 1;
cfg.mult_comp_corr = 1;

% Smoothing
cfg.doSmooth               = 1;
cfg.spine_smooth_fwhm_mm   = 20;
cfg.brain_smooth_fwhm_mm   = 8;
cfg.spine_smooth_radius_mm = 3 * (cfg.spine_smooth_fwhm_mm  / 2.355);
cfg.brain_smooth_radius_mm = 3 * (cfg.brain_smooth_fwhm_mm / 2.355);

% Rectified results file (from main pipeline) for spine comparison plot
cfg.rect_results_file = fullfile(cfg.save_dir, 'subResult_subOP00212_BS.mat');
cfg.p1_idx   = 1;   % index of OP00212 in subjResults
cfg.saveFigs = 1;

% Steps to run
run_step1 = 1;   % Brain-EMG
run_step2 = 1;   % Spine-EMG

%% =========================================================================
%  SETUP
%% =========================================================================

addpath(cfg.bsc_source_path)
addpath(cfg.spm_path)
spm('defaults','EEG')
addpath(cfg.fieldtrip_path)
ft_defaults

rng(cfg.rng_seed)

if cfg.saveFigs && ~exist(cfg.fig_dir,'dir'), mkdir(cfg.fig_dir); end

fprintf('\n=== UNRECTIFIED EMG SENSITIVITY - %s ===\n', cfg.sub)
fprintf('EMG rectification: OFF (unrectified)\n\n')

%% =========================================================================
%  LOAD DATA (shared across both steps)
%% =========================================================================

run      = '001';
datafile = fullfile(cfg.data_root, ['sub-' cfg.sub], 'ses-001', 'meg', ...
    sprintf('pmergedoe1000mspddfflo45hi45hfcstatic_%s_array1.mat', run));

D       = spm_eeg_load(datafile);
grad_mm = D.sensors('MEG');
ftdat   = spm2fieldtrip(D);

% Remove bad channels
badchans = D.chanlabels(D.badchannels);
cfg_ft   = []; cfg_ft.channel = setdiff(ftdat.label, badchans);
ftdat    = ft_selectdata(cfg_ft, ftdat);

% No EMG rectification — ftdat used as-is

cmap = make_cmap(cfg.fieldtrip_path);

%% =========================================================================
%  STEP 1 — Brain-EMG coherence (unrectified EMG)
%  Uses all channels for leadfield (grad_mm, not braingrad)
%% =========================================================================

if run_step1
    fprintf('>>> STEP 1: Brain-EMG (unrectified)\n')

    %% Load brain geometry
    t          = load(cfg.geomfile);
    mesh_brain = t.mesh_brain;
    mesh_brain.unit = 'mm';

    %% Head model
    cfg_ft = []; cfg_ft.method = 'singleshell';
    vol = ft_prepare_headmodel(cfg_ft, mesh_brain);

    %% Leadfield — all channels (grad_mm)
    cfg_ft            = [];
    cfg_ft.sourcemodel = t.sources_brain;
    cfg_ft.headmodel   = vol;
    cfg_ft.grad        = grad_mm;   % all channels, not braingrad
    cfg_ft.reducerank  = 'no';
    LF_brain = ft_prepare_leadfield(cfg_ft);


    %% Frequency data (full ftdat — channel matching handled by FieldTrip)
    [statdat_br, restdat_br, freqdat_br] = make_freq_data(ftdat, ftdat.trialinfo, cfg.fband);

   %% Common spatial filter
cfg_ft                        = [];
cfg_ft.grid                   = t.sources_brain;
cfg_ft.headmodel              = vol;
cfg_ft.sourcemodel.leadfield  = LF_brain.leadfield;
cfg_ft.dics.keepfilter        = 'yes';
cfg_ft.dics.lambda            = '10%';
cfg_ft.method                 = 'dics';
cfg_ft.refchan                = 'EXG1';
coh_source_br = ft_sourceanalysis(cfg_ft, freqdat_br);

%% Permutation test
cfg_ft                        = [];
cfg_ft.grid                   = t.sources_brain;
cfg_ft.headmodel              = vol;
cfg_ft.sourcemodel.leadfield  = LF_brain.leadfield;
cfg_ft.dics.filter            = coh_source_br.avg.filter;
cfg_ft.dics.lambda            = '10%';
cfg_ft.method                 = 'dics';
cfg_ft.refchan                = 'EXG1';
cfg_ft.permutation            = 'yes';
cfg_ft.numpermutation         = cfg.numpermutation;
source_perm_br = ft_sourceanalysis(cfg_ft, statdat_br, restdat_br);

    nsourcepoints_br = size(source_perm_br.pos, 1);
    nPerm_br         = numel(source_perm_br.trialA);
    [coh_diff_br, cohDiff_perm_br] = extract_coh_diff(source_perm_br, nsourcepoints_br, nPerm_br);

    %% Smoothing (inside sources only)
    if cfg.doSmooth
        inside_idx_br   = find(source_perm_br.inside);
        inside_pos_br   = source_perm_br.pos(inside_idx_br, :);
        Wsm_br          = make_gaussian_smoother(inside_pos_br, cfg.brain_smooth_fwhm_mm, cfg.brain_smooth_radius_mm);
        cd_in           = coh_diff_br(inside_idx_br);
        cp_in           = cohDiff_perm_br(inside_idx_br, :);
        cd_in           = Wsm_br * cd_in;
        cp_in           = Wsm_br * cp_in;
        coh_diff_br(inside_idx_br)       = cd_in;
        cohDiff_perm_br(inside_idx_br,:) = cp_in;
    end

    %% Threshold and stats
    thr95_br    = compute_threshold(cohDiff_perm_br(inside_idx_br,:), cfg.mult_comp_corr);
    pvals_br    = compute_pvals(coh_diff_br, cohDiff_perm_br);
    invp_br     = -log10(pvals_br);
    mask_br     = coh_diff_br > thr95_br;
    invp_masked_br = invp_br; invp_masked_br(~mask_br) = 0;
    invpthr     = -log10(0.05);

    fprintf('  max coh diff = %.4f, thr95 = %.4f, n_sig = %d\n', ...
        max(coh_diff_br), thr95_br, sum(mask_br))

    %% MNI peak
    [maxpos_br, x_mni_br] = get_mni_max(coh_diff_br, source_perm_br, cfg.T_mat_path);
    fprintf('  Peak MNI: [%.1f %.1f %.1f]\n', x_mni_br(1,:))

    %% Interpolate and plot
    [brain_int_br, source_mask_int_br] = interpolate_brain(coh_source_br, invp_masked_br, mask_br, mesh_brain);
    hfig_br = plot_brain_surface(brain_int_br, mask_br, invp_br, invpthr, cmap, mesh_brain, ...
        sprintf('%s brain-EMG coherence (unrectified)', cfg.sub), ...
        cfg.doSmooth, cfg.brain_smooth_fwhm_mm);

    %% Save
    subjResults_brain(1).subjID     = cfg.sub;
    subjResults_brain(1).coh_diff   = coh_diff_br;
    subjResults_brain(1).sig_mask   = double(mask_br);
    subjResults_brain(1).pos        = source_perm_br.pos;
    subjResults_brain(1).inside     = source_perm_br.inside;
    subjResults_brain(1).thr95      = thr95_br;
    subjResults_brain(1).maxpos     = maxpos_br;
    subjResults_brain(1).maxpos_mni = x_mni_br;
    subjResults_brain(1).brain_int  = brain_int_br;
    save(fullfile(cfg.save_dir, sprintf('brainResult_sub%s_unrectified.mat', cfg.sub)), ...
        'subjResults_brain')
    fprintf('  Brain results saved.\n')

    if cfg.saveFigs
        savefig(hfig_br, fullfile(cfg.fig_dir, ...
            sprintf('step1_brainEMG_unrectified_%s.fig', cfg.sub)));
        saveas(hfig_br,  fullfile(cfg.fig_dir, ...
            sprintf('step1_brainEMG_unrectified_%s.png', cfg.sub)));
    end
    fprintf('>>> STEP 1 complete.\n\n')
end

%% =========================================================================
%  STEP 2 — Spine-EMG coherence (unrectified EMG, BSLaw forward model)
%% =========================================================================

if run_step2
    fprintf('>>> STEP 2: Spine-EMG (unrectified)\n')

    %% Load spine geometry
    sources_cent = t.sources_cent;
    mesh_torso   = t.mesh_torso;
    mesh_wm      = t.mesh_wm;
    mesh_wm.unit = 'mm';

    nsourcepoints = size(sources_cent.pos, 1);

    %% Load BSLaw leadfield and match channels
    lf_data = load(cfg.lf_bs_path);
    lf_raw  = lf_data.leadfield_bs;

    data_meg_labels        = ftdat.label(~strcmp(ftdat.label,'EXG1'));
    [common_labels, idx_lf] = intersect(lf_raw.label, data_meg_labels, 'stable');
    fprintf('  Data MEG: %d  |  LF: %d  |  Matched: %d\n', ...
        numel(data_meg_labels), numel(lf_raw.label), numel(common_labels));
    if numel(common_labels) < numel(data_meg_labels)
        fprintf('  WARNING: %d MEG channels not in leadfield\n', ...
            numel(data_meg_labels) - numel(common_labels));
    end

    Lf        = lf_raw;
    Lf.label  = common_labels;
    Lf.pos    = sources_cent.pos;
    Lf.inside = ones(nsourcepoints, 1);
    for i = 1:numel(lf_raw.leadfield)
        if ~isempty(lf_raw.leadfield{i})
            Lf.leadfield{i} = lf_raw.leadfield{i}(idx_lf, :);
        end
    end

    %% Volume conductor (infinite homogeneous, matches main spine script)
    cfg_ft            = [];
    cfg_ft.method     = 'infinite';
    cfg_ft.siunits    = 1;
    cfg_ft.grad       = grad_mm;
    cfg_ft.conductivity = 1;
    dummyvol = ft_prepare_headmodel(cfg_ft, mesh_torso);

    %% Frequency data — restrict to leadfield channels + EXG1
    [statdat_sp, restdat_sp, freqdat_sp] = make_freq_data(ftdat, ftdat.trialinfo, cfg.fband);

    cfg_sel          = [];
    cfg_sel.channel  = [Lf.label; {'EXG1'}];   % CRITICAL: must match leadfield
    statdat_sp  = ft_selectdata(cfg_sel, statdat_sp);
    restdat_sp  = ft_selectdata(cfg_sel, restdat_sp);
    freqdat_sp  = ft_selectdata(cfg_sel, freqdat_sp);

    %% Sourcemodel struct
    sourcemodel_sp            = [];
    sourcemodel_sp.pos        = Lf.pos;
    sourcemodel_sp.unit       = 'mm';
    sourcemodel_sp.inside     = logical(Lf.inside);
    sourcemodel_sp.leadfield  = Lf.leadfield;
    sourcemodel_sp.label      = Lf.label;

    %% Common spatial filter
    cfg_ft                    = [];
    cfg_ft.sourcemodel        = sourcemodel_sp;
    cfg_ft.headmodel          = dummyvol;
    cfg_ft.dics.keepfilter    = 'yes';
    cfg_ft.dics.lambda        = '10%';
    cfg_ft.method             = 'dics';
    cfg_ft.refchan            = 'EXG1';
    coh_source_sp = ft_sourceanalysis(cfg_ft, freqdat_sp);

    %% Permutation test
    cfg_ft                    = [];
    cfg_ft.sourcemodel        = sourcemodel_sp;
    cfg_ft.headmodel          = dummyvol;
    cfg_ft.dics.filter        = coh_source_sp.avg.filter;
    cfg_ft.dics.lambda        = '10%';
    cfg_ft.method             = 'dics';
    cfg_ft.refchan            = 'EXG1';
    cfg_ft.permutation        = 'yes';
    cfg_ft.numpermutation     = cfg.numpermutation;
    source_perm_sp = ft_sourceanalysis(cfg_ft, statdat_sp, restdat_sp);

    nPerm_sp = numel(source_perm_sp.trialA);
    [coh_diff_sp, cohDiff_perm_sp] = extract_coh_diff(source_perm_sp, nsourcepoints, nPerm_sp);

    %% Smoothing
    if cfg.doSmooth
        Wsm_sp       = make_gaussian_smoother(sources_cent.pos, ...
            cfg.spine_smooth_fwhm_mm, cfg.spine_smooth_radius_mm);
        cohDiff_perm_sp = Wsm_sp * cohDiff_perm_sp;
        coh_diff_sp     = Wsm_sp * coh_diff_sp;
    end

    %% Threshold and stats
    thr95_sp    = compute_threshold(cohDiff_perm_sp, cfg.mult_comp_corr);
    pvals_sp    = compute_pvals(coh_diff_sp, cohDiff_perm_sp);
    invp_sp     = -log10(pvals_sp);
    mask_sp     = coh_diff_sp > thr95_sp;
    invp_smooth_sp = smooth_invp(coh_diff_sp, cohDiff_perm_sp, nsourcepoints, nPerm_sp);

    fprintf('  Unrectified: max coh diff = %.4f, thr95 = %.4f, n_sig = %d\n', ...
        max(coh_diff_sp), thr95_sp, sum(mask_sp))

    %% Null maxima histogram (P1 only)
    [~, maxIdx_perm_sp] = max(cohDiff_perm_sp, [], 1);
    [~, obsMaxIdx_sp]   = max(coh_diff_sp);
    cord_pos            = sources_cent.pos(:, 2);

    hfig_null = figure('Color','w','Position',[100 100 600 450]);
    hold on;
    histogram(cord_pos(maxIdx_perm_sp), 44, ...
        'FaceColor',[0.75 0.75 0.75],'EdgeColor','k','LineWidth',0.8);
    xline(cord_pos(obsMaxIdx_sp), '-', 'Color',[0.2 0 0], 'LineWidth', 2);
    xlabel('Cranio-caudal position (mm)','FontSize',14);
    ylabel('Count','FontSize',14);
    legend({'Null maxima','Observed maximum'},'Location','best','FontSize',14,'Box','off');
    set(gca,'FontSize',14,'LineWidth',1.2,'TickDir','out'); box off;
    title(sprintf('%s — null maxima spine-EMG (unrectified, BS)', cfg.sub),'Interpreter','none');
    if cfg.saveFigs
        savefig(hfig_null, fullfile(cfg.fig_dir, ...
            sprintf('step2_null_maxima_unrectified_%s.fig', cfg.sub)));
    end
    close(hfig_null);

    %% Save spine results
    [~, obsIdx_sp] = max(invp_smooth_sp);
    subjResults_spine(1).subjID      = cfg.sub;
    subjResults_spine(1).coh_diff    = coh_diff_sp;
    subjResults_spine(1).thr95       = thr95_sp;
    subjResults_spine(1).sig_mask    = mask_sp;
    subjResults_spine(1).pos         = sources_cent.pos;
    subjResults_spine(1).inside      = sources_cent.inside;
    subjResults_spine(1).maxdiff.idx = obsIdx_sp;
    subjResults_spine(1).maxdiff.pos = sources_cent.pos(obsIdx_sp, :);
    subjResults_spine(1).invp_smooth = invp_smooth_sp;
    save(fullfile(cfg.save_dir, sprintf('spineResult_sub%s_unrectified.mat', cfg.sub)), ...
        'subjResults_spine')
    fprintf('  Spine results saved.\n')

    %% Load rectified result for comparison
    rect_coh_diff = [];
    rect_thr95    = NaN;
    rect_mask     = false(nsourcepoints, 1);
    if exist(cfg.rect_results_file, 'file')
        rect_res  = load(cfg.rect_results_file, 'subjResults');
        sr        = rect_res.subjResults(cfg.p1_idx);
        rect_coh_diff = sr.coh_diff;
        rect_thr95    = sr.thr95;
        rect_mask     = logical(sr.sig_mask);
        fprintf('  Rectified:   max coh diff = %.4f, thr95 = %.4f, n_sig = %d\n', ...
            max(rect_coh_diff), rect_thr95, sum(rect_mask))
    else
        warning('Rectified results file not found: %s', cfg.rect_results_file)
    end

    %% Comparison line plot — rectified vs unrectified
    col_rect   = [55  126 184] / 255;   % blue
    col_unrect = [228  26  28] / 255;   % red
    x          = sources_cent.pos(:, 2);

    hfig_sp = figure('Color','w','Position',[100 100 560 500]);
    hold on;

    % Rectified trace
    if ~isempty(rect_coh_diff)
        for i = 1:length(x)-1
            if rect_mask(i) && rect_mask(i+1)
                plot(x(i:i+1), rect_coh_diff(i:i+1), '-', ...
                    'Color', col_rect, 'LineWidth', 1.5, 'HandleVisibility','off')
            else
                plot(x(i:i+1), rect_coh_diff(i:i+1), '-', ...
                    'Color', [0.75 0.75 0.75], 'LineWidth', 0.8, 'HandleVisibility','off')
            end
        end
        plot(x(rect_mask), rect_coh_diff(rect_mask), '.', ...
            'Color', col_rect, 'MarkerSize', 12, 'HandleVisibility','off')
        h_rect = plot(nan, nan, '-', 'Color', col_rect, 'LineWidth', 2, ...
            'DisplayName','Rectified');
        xline(rect_thr95, '--', 'Color', col_rect, 'LineWidth', 0.8, ...
            'HandleVisibility','off');
    end

    % Unrectified trace
    for i = 1:length(x)-1
        if mask_sp(i) && mask_sp(i+1)
            plot(x(i:i+1), coh_diff_sp(i:i+1), '-', ...
                'Color', col_unrect, 'LineWidth', 1.5, 'HandleVisibility','off')
        else
            plot(x(i:i+1), coh_diff_sp(i:i+1), '-', ...
                'Color', [0.85 0.75 0.75], 'LineWidth', 0.8, 'HandleVisibility','off')
        end
    end
    plot(x(mask_sp), coh_diff_sp(mask_sp), '.', ...
        'Color', col_unrect, 'MarkerSize', 12, 'HandleVisibility','off')
    h_unrect = plot(nan, nan, '-', 'Color', col_unrect, 'LineWidth', 2, ...
        'DisplayName','Unrectified');
    xline(thr95_sp, '--', 'Color', col_unrect, 'LineWidth', 0.8, ...
        'HandleVisibility','off');

    yline(0, ':k', 'HandleVisibility','off')

    if ~isempty(rect_coh_diff)
        legend([h_rect h_unrect], 'Location','best','Box','off','FontSize',12);
    else
        legend(h_unrect, 'Location','best','Box','off','FontSize',12);
    end

    xlabel('Cranial-caudal position (mm)', 'FontSize', 13);
    ylabel('Coherence difference (contraction - rest)', 'FontSize', 13);
    title(sprintf('%s spine-EMG: rectified vs unrectified (BSLaw)', cfg.sub), ...
        'Interpreter','none','FontSize',13);
    set(gca,'FontSize',13); box off;

    if cfg.saveFigs
        savefig(hfig_sp, fullfile(cfg.fig_dir, ...
            sprintf('step2_spineEMG_rect_vs_unrect_%s.fig', cfg.sub)));
        print(hfig_sp, fullfile(cfg.fig_dir, ...
            sprintf('step2_spineEMG_rect_vs_unrect_%s.png', cfg.sub)), '-dpng', '-r300');
    end
    fprintf('>>> STEP 2 complete.\n\n')
end

fprintf('\n=== UNRECTIFIED EMG ANALYSIS FINISHED ===\n\n')


%% =========================================================================
%  LOCAL UTILITY FUNCTIONS
%% =========================================================================

function [statdat, restdat, freqdat] = make_freq_data(dat, trialinfo, fband)
% Compute CSD freq data with balanced trials for stat/rest and combined.
    cfg = []; cfg.output = 'powandcsd'; cfg.method = 'mtmfft';
    cfg.foilim = fband; cfg.tapsmofrq = 1; cfg.keeptrials = 'yes';
    freqdat_tr = ft_freqanalysis(cfg, dat);
    cfg = []; cfg.avgoverfreq = 'yes';
    freqdat_tr = ft_selectdata(cfg, freqdat_tr);

    statidx = find(trialinfo == 1);
    restidx = find(trialinfo == 2);
    nTrials = min(numel(statidx), numel(restidx));

    cfg = []; cfg.trials = statidx(1:nTrials);
    statdat = ft_selectdata(cfg, freqdat_tr);
    cfg = []; cfg.trials = restidx(1:nTrials);
    restdat = ft_selectdata(cfg, freqdat_tr);

    cfg = []; cfg.output = 'powandcsd'; cfg.method = 'mtmfft';
    cfg.foilim = fband; cfg.tapsmofrq = 1; cfg.keeptrials = 'no';
    freqdat = ft_freqanalysis(cfg, dat);
    cfg = []; cfg.avgoverfreq = 'yes';
    freqdat = ft_selectdata(cfg, freqdat);
end

% -------------------------------------------------------------------------
function [coh_diff, cohDiff_perm] = extract_coh_diff(source_perm, nsourcepoints, nPerm)
    cohDiff_perm = zeros(nsourcepoints, nPerm);
    for i = 1:nPerm
        cohDiff_perm(:,i) = source_perm.trialA(i).coh - source_perm.trialB(i).coh;
    end
    coh_diff = source_perm.avgA.coh - source_perm.avgB.coh;
end

% -------------------------------------------------------------------------
function thr = compute_threshold(cohDiff_perm, mult_comp_corr)
    if mult_comp_corr
        thr = prctile(max(cohDiff_perm, [], 1), 95);
    else
        thr = prctile(cohDiff_perm, 95, 2);
        warning('Uncorrected per-source threshold used.');
    end
end

% -------------------------------------------------------------------------
function pvals = compute_pvals(coh_diff, cohDiff_perm)
    nPerm = size(cohDiff_perm, 2);
    pvals = (sum(cohDiff_perm >= coh_diff, 2) + 1) / (nPerm + 1);
end

% -------------------------------------------------------------------------
function [brain_int, source_mask_int] = interpolate_brain(coh_source, invp_masked, mask, mesh_brain)
    source_p         = coh_source;
    source_p.avg.coh = invp_masked;
    cfg = []; cfg.parameter = 'coh'; cfg.interpmethod = 'sphere_avg';
    cfg.sphereradius = 10;
    brain_int = ft_sourceinterpolate(cfg, source_p, mesh_brain);
    source_mask         = coh_source;
    source_mask.avg.coh = double(mask);
    source_mask_int = ft_sourceinterpolate(cfg, source_mask, mesh_brain);
    brain_int.mask  = source_mask_int.coh;
end

% -------------------------------------------------------------------------
function hfig = plot_brain_surface(brain_int, mask, invp, invpthr, cmap, mesh_brain, ttl, doSmooth, fwhm)
    hfig = figure;
    cfg = []; cfg.figure = 'gcf'; cfg.method = 'surface';
    cfg.funparameter = 'coh'; cfg.funcolormap = cmap;
    if any(mask)
        cfg.funcolorlim = [invpthr max(invp(mask))];
    else
        cfg.funcolorlim = [invpthr invpthr+1];
    end
    cfg.projmethod = 'nearest'; cfg.surffile = mesh_brain;
    ft_sourceplot(cfg, brain_int);
    view(176,-10); camlight; ax = gca; ax.FontSize = 14;
    hpatch = findobj(gcf,'Type','patch'); set(hpatch,'FaceAlpha',0.9)
    if doSmooth
        title(sprintf('%s (smoothed %d mm)', ttl, fwhm),'Interpreter','none')
    else
        title(ttl,'Interpreter','none')
    end
end

% -------------------------------------------------------------------------
function cmap = make_cmap(fieldtrip_path)
    ncol = 256;
    addpath(fullfile(fieldtrip_path,'external','matplotlib'))
    cmap = [[0.92 0.92 0.92]; flipud(magma(ncol-1))];
end

% -------------------------------------------------------------------------
function W = make_gaussian_smoother(pos_mm, fwhm_mm, radius_mm)
    sigma = fwhm_mm / 2.355;
    if nargin < 3 || isempty(radius_mm), radius_mm = 3*sigma; end
    N   = size(pos_mm, 1);
    Mdl = KDTreeSearcher(pos_mm);
    [idx, dist] = rangesearch(Mdl, pos_mm, radius_mm);
    ii = []; jj = []; vv = [];
    for i = 1:N
        j = idx{i}; d = dist{i};
        w = exp(-0.5*(d./sigma).^2);
        ii = [ii; repmat(i,numel(j),1)]; jj = [jj; j(:)]; vv = [vv; w(:)];
    end
    W  = sparse(ii,jj,vv,N,N);
    rs = full(sum(W,2)); rs(rs==0) = 1;
    W  = spdiags(1./rs,0,N,N) * W;
end

% -------------------------------------------------------------------------
function [maxpos, x_mni] = get_mni_max(coh_diff, source_perm, T_mat_path)
    maxval   = max(coh_diff);
    maxidx   = find(coh_diff == maxval);
    maxpos   = source_perm.pos(maxidx, :);
    load(T_mat_path, 'T'); T_inv = inv(T);
    maxpos_h = [maxpos, ones(size(maxpos,1),1)]';
    x_mni    = (T_inv * maxpos_h)'; x_mni = x_mni(:,1:3);
    if length(maxidx) > 1, disp('  Multiple max locations'); end
end

% -------------------------------------------------------------------------
function invp_smooth = smooth_invp(coh_diff, cohDiff_perm, nsourcepoints, nPerm)
    invp_smooth = zeros(nsourcepoints, 1);
    for s = 1:nsourcepoints
        permDist = sort(cohDiff_perm(s,:));
        obsVal   = coh_diff(s);
        xgrid    = linspace(min(permDist), max(permDist), 200);
        p_emp    = arrayfun(@(x) (sum(permDist >= x)+1)/(nPerm+1), xgrid);
        logp_sm  = smooth(xgrid, -log10(p_emp), 0.15, 'loess');
        obsVal_c = min(max(obsVal, xgrid(1)), xgrid(end));
        invp_smooth(s) = interp1(xgrid, logp_sm, obsVal_c, 'linear');
    end
end

% -------------------------------------------------------------------------
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
