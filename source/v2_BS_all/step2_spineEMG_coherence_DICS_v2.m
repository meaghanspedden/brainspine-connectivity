%% step2_spineEMG_coherence_DICS_v2.m
% Step 2 (v2): Spine-EMG DICS coherence with Biot-Savart (BS) forward model.
%
% Uses the BS leadfield with label-based channel matching. Lambda = 10%,
% smoothing 20mm FWHM. All 9 spine subjects. Outputs saved with _BS
% suffix to distinguish from older (BEM) results.
%
% Differs from step2_spineEMG_coherence_DICS.m (v1) in its orientation
% analysis: v1 produces a per-subject dot-plot of observed vs. null
% orientation (median +/- IQR). This version (v2) instead produces:
%   - Per-subject null maxima location histogram
%   - Dominant orientation extraction per subject (SVD at peak source)
%   - Participant-1-only null distribution stratified by orientation axis
%   - Group-level dominant orientation bar chart across participants ("Figure C")
% Run both and compare outputs to decide which orientation analysis to keep.
%
% Inputs (must exist before running):
%   geomfile       - geometries_experimental.mat: sources_cent, mesh_torso,
%                     mesh_wm, mesh_bone, mesh_lungs, mesh_heart
%   geomfile_brain - geometries_cervical_realistic.mat: mesh_brain only
%                     (TODO: mesh_brain still needs adding to geomfile
%                     directly; see warning below)
%   lf_path        - leadfield_experimental_bslaw_experimental.mat: leadfield_bs
%   data_root      - per-subject SPM .mat files (see datafile pattern below)
%
% Outputs (written to save_dir / fig_dir):
%   subResult_sub<ID>_BS.mat                              - per-subject coh_diff, thr95, mask, invp_smooth, pvals
%   groupRes_spine_DICS__BS.mat                            - subjResults struct for all subjects (incl. dom_ori)
%   cluster_spineEMG_pos_BS.mat                            - VE ROI cluster positions
%   step2_sub<ID>_spineEMG_coherence_BS.fig                - per-subject coherence map (with torso)
%   step2_sub1_spineEMG_meshonly_BS.fig                    - Participant 1 mesh-only map
%   step2_sub<ID>_null_maxima_BS.fig                       - per-subject null maxima location histogram
%   step2_subOP00212_null_distribution_by_orientation_BS.fig - P1-only null dist. by orientation axis
%   step2_group_spineEMG_prevalence_BS.fig                 - group prevalence map (with torso)
%   step2_group_spineEMG_prevalence_meshonly_BS.fig        - group prevalence, mesh only
%   step2_group_spineEMG_subject_lines_BS.fig              - per-subject coherence-difference traces
%   step2_group_dominant_orientation_BS.fig / .png         - Figure C: group orientation bar chart
%   step2_spinal_prevalence_VE_cluster_BS.fig              - prevalence + VE cluster ROI
%
% rng(1) is set below to fix the permutation test seed for reproducibility.

clear all; close all; clc;

%% =========================================================================
%  USER CONFIG
%% =========================================================================
fieldtrip_path = 'C:\Users\mspedden\Documents\fieldtrip';
spm_path       = 'C:\Users\mspedden\Documents\spm';
bsc_path       = 'C:\Users\mspedden\Documents\brainspineconnectivity\source';
data_root      = 'C:\spinecoh_data';
save_dir       = 'C:\Users\mspedden\Documents\brainspine_savetest';

geomfile       = 'C:\Leadfields meshes\geometries_experimental.mat';
geomfile_brain = 'C:\Leadfields meshes\geometries_cervical_realistic.mat';
lf_path        = 'C:\Leadfields meshes\leadfield_experimental_bslaw_experimental.mat';

subs_spine = {'OP00212','OP00213','OP00215','OP00219', ...
              'OP00220','OP00221','OP00224','OP00225','OP00226'};

fband          = [10 35];
numpermutation = 500;
mult_comp_corr = 1;
lambda         = '10%';
fwhm_mm        = 20;
radius_mm      = 3 * (fwhm_mm / 2.355);
out_suffix     = '_BS';
rng(1);   % fix seed for reproducible permutation test

%% =========================================================================
%  SETUP
%% =========================================================================
addpath(bsc_path); addpath(spm_path);
spm('defaults','EEG');
addpath(fieldtrip_path);
ft_defaults;

if ~exist(save_dir,'dir'), mkdir(save_dir); end
fig_dir = fullfile(save_dir, 'figures');
if ~exist(fig_dir,'dir'), mkdir(fig_dir); end

%% =========================================================================
%  LOAD GEOMETRY
%% =========================================================================
fprintf('Loading geometry...\n');
geom_exp = load(geomfile);

sources_cent = geom_exp.sources_cent;
mesh_torso   = geom_exp.mesh_torso;
mesh_wm      = geom_exp.mesh_wm;
mesh_bone    = geom_exp.mesh_bone;
mesh_lungs   = geom_exp.mesh_lungs;
mesh_heart   = geom_exp.mesh_heart;

warning('Need to add brain mesh to new geom!!!!')
geom_brain = load(geomfile_brain, 'mesh_brain');
mesh_brain = geom_brain.mesh_brain;

mesh_wm.unit = 'mm';

nsourcepoints = size(sources_cent.pos, 1);
fprintf('  Source space: %d points\n', nsourcepoints);

%% =========================================================================
%  BUILD SMOOTHER ONCE
%% =========================================================================
fprintf('Building Gaussian smoother (FWHM=%d mm)...\n', fwhm_mm);
Wsm = make_gaussian_smoother(sources_cent.pos, fwhm_mm, radius_mm);
nnz_per_row = full(sum(Wsm > 0, 2));
selfw       = full(diag(Wsm));
fprintf('  Neighbours/row: median %.1f (min %d, max %d)\n', ...
    median(nnz_per_row), min(nnz_per_row), max(nnz_per_row));
fprintf('  Self-weight:    median %.3f (min %.3f, max %.3f)\n\n', ...
    median(selfw), min(selfw), max(selfw));

%% =========================================================================
%  COLOURMAP
%% =========================================================================
ncol     = 256;
addpath(fullfile(fieldtrip_path,'external','matplotlib'));
cmap_hot = [[0.92 0.92 0.92]; flipud(magma(ncol-1))];

%% =========================================================================
%  SUBJECT LOOP
%% =========================================================================
subjResults = struct();

for ss = 1:length(subs_spine)
    sub = subs_spine{ss};
    fprintf('\n========================================\n');
    fprintf('  Subject %s (%d/%d)\n', sub, ss, length(subs_spine));
    fprintf('========================================\n');
    t_sub = tic;

    %% Load data
    run = '001'; if strcmp(sub,'OP00224'), run = '002'; end
    datafile = fullfile(data_root, ['sub-' sub], 'ses-001', 'meg', ...
        sprintf('pmergedoe1000mspddfflo45hi45hfcstatic_%s_array1.mat', run));

    D       = spm_eeg_load(datafile);
    grad_mm = D.sensors('MEG');
    ftdat   = spm2fieldtrip(D);

    % Remove bad channels
    badchans = D.chanlabels(D.badchannels);
    cfg = []; cfg.channel = setdiff(ftdat.label, badchans);
    ftdat = ft_selectdata(cfg, ftdat);

    % Rectify EMG
    cfg = []; cfg.rectify = 'yes'; cfg.channel = 'EXG1';
    ftdatr = ft_preprocessing(cfg, ftdat);
    for k = 1:length(ftdat.trial)
        ftdat.trial{k}(end,:) = ftdatr.trial{k};
    end

    %% Load and build leadfield — label-based matching
    fprintf('  Building leadfield...\n');
    lf_data = load(lf_path);
    lf_raw  = lf_data.leadfield_bs;

    data_meg_labels        = ftdat.label(~strcmp(ftdat.label,'EXG1'));
    [common_labels,idx_lf] = intersect(lf_raw.label, data_meg_labels, 'stable');
    fprintf('  Data MEG: %d  |  LF: %d  |  Matched: %d\n', ...
        numel(data_meg_labels), numel(lf_raw.label), numel(common_labels));
    if numel(common_labels) < numel(data_meg_labels)
        fprintf('  WARNING: %d channels not in leadfield\n', ...
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

    %% Volume conductor
    cfg = []; cfg.method = 'infinite'; cfg.siunits = 1;
    cfg.grad = grad_mm; cfg.conductivity = 1;
    dummyvol = ft_prepare_headmodel(cfg, mesh_torso);

    %% Frequency data
    cfg_fr  = []; cfg_fr.output = 'powandcsd'; cfg_fr.method = 'mtmfft';
    cfg_fr.foilim = fband; cfg_fr.tapsmofrq = 1; cfg_fr.keeptrials = 'yes';
    cfg_av  = []; cfg_av.avgoverfreq = 'yes';
    cfg_sel = []; cfg_sel.channel = [Lf.label; {'EXG1'}];

    freqdat_tr = ft_freqanalysis(cfg_fr, ftdat);
    freqdat_tr = ft_selectdata(cfg_av, freqdat_tr);

    trialinfo = ftdat.trialinfo;
    statidx   = find(trialinfo == 1);
    restidx   = find(trialinfo == 2);
    nTrials   = min(numel(statidx), numel(restidx));
    fprintf('  nTrials (stat/rest): %d\n', nTrials);

    cfg = []; cfg.trials = statidx(1:nTrials);
    statdat = ft_selectdata(cfg, freqdat_tr);
    cfg = []; cfg.trials = restidx(1:nTrials);
    restdat = ft_selectdata(cfg, freqdat_tr);

    cfg_fr2 = []; cfg_fr2.output = 'powandcsd'; cfg_fr2.method = 'mtmfft';
    cfg_fr2.foilim = fband; cfg_fr2.tapsmofrq = 1; cfg_fr2.keeptrials = 'no';
    freqdat = ft_freqanalysis(cfg_fr2, ftdat);
    cfg_av2 = []; cfg_av2.avgoverfreq = 'yes';
    freqdat = ft_selectdata(cfg_av2, freqdat);

    statdat = ft_selectdata(cfg_sel, statdat);
    restdat = ft_selectdata(cfg_sel, restdat);
    freqdat = ft_selectdata(cfg_sel, freqdat);

    %% Sourcemodel
    sourcemodel           = [];
    sourcemodel.pos       = Lf.pos;
    sourcemodel.unit      = 'mm';
    sourcemodel.inside    = logical(Lf.inside);
    sourcemodel.leadfield = Lf.leadfield;
    sourcemodel.label     = Lf.label;

    %% Common spatial filter
    fprintf('  Computing common spatial filter...\n');
    cfg_dics                 = [];
    cfg_dics.sourcemodel     = sourcemodel;
    cfg_dics.headmodel       = dummyvol;
    cfg_dics.dics.keepfilter = 'yes';
    cfg_dics.dics.lambda     = lambda;
    cfg_dics.method          = 'dics';
    cfg_dics.refchan         = 'EXG1';
    coh_source = ft_sourceanalysis(cfg_dics, freqdat);

    %% Permutation test
    fprintf('  Running permutation test (%d permutations)...\n', numpermutation);
    cfg_perm                    = [];
    cfg_perm.sourcemodel        = sourcemodel;
    cfg_perm.headmodel          = dummyvol;
    cfg_perm.dics.filter        = coh_source.avg.filter;
    cfg_perm.dics.lambda        = lambda;
    cfg_perm.method             = 'dics';
    cfg_perm.refchan            = 'EXG1';
    cfg_perm.permutation        = 'yes';
    cfg_perm.numpermutation     = numpermutation;
    source_perm = ft_sourceanalysis(cfg_perm, statdat, restdat);

    nPerm = numel(source_perm.trialA);
    [coh_diff, cohDiff_perm] = extract_coh_diff(source_perm, nsourcepoints, nPerm);

    %% Smoothing
    cohDiff_perm = Wsm * cohDiff_perm;
    coh_diff     = Wsm * coh_diff;

    %% Threshold
    thr95 = compute_threshold(cohDiff_perm, mult_comp_corr, nsourcepoints);
    mask  = coh_diff > thr95;

    fprintf('  Threshold (FWE p<0.05): %.6f\n', thr95);
    fprintf('  Significant sources:    %d / %d\n', sum(mask), nsourcepoints);
    [peak_val, peak_idx] = max(coh_diff);
    fprintf('  Peak coh_diff: %.4e at y=%.1f mm (source %d)\n', ...
        peak_val, sources_cent.pos(peak_idx,2), peak_idx);
    fprintf('  Subject time: %.1f min\n', toc(t_sub)/60);

    %% p-values and smoothed inverse p
    pvals       = compute_pvals(coh_diff, cohDiff_perm);
    invp        = -log10(pvals);
    invp_masked = invp; invp_masked(~mask) = 0;
    invpthr     = -log10(0.05);
    invp_smooth = smooth_invp(coh_diff, cohDiff_perm, nsourcepoints, nPerm);

    %% Null maxima location and max per permutation (needed for orientation figures)
    maxPerm          = max(cohDiff_perm, [], 1);   % [1 x nPerm] — max over sources
    [~, maxIdx_perm] = max(cohDiff_perm, [], 1);   % index of max source per perm
    xpos             = sources_cent.pos(:,2);
    [~, obsMaxIdx]   = max(coh_diff);

    %% -----------------------------------------------------------------------
    %  DOMINANT ORIENTATION AT PEAK SOURCE (all subjects)
    %  SVD of DICS filter at peak coh_diff source -> leading left singular
    %  vector -> decompose into x/y/z components
    %% -----------------------------------------------------------------------
    filters = coh_source.avg.filter;   % cell {nsources x 1}, each [3 x nchans]
    F_peak  = filters{peak_idx};
    if ~isempty(F_peak)
        [U,~,~] = svd(F_peak, 'econ');
        dom_vec  = abs(U(:,1));         % |components| of dominant orientation
        subjResults(ss).dom_ori = dom_vec;   % [3x1]: [x; y; z]
    else
        subjResults(ss).dom_ori = [nan; nan; nan];
    end

    %% Null maxima location histogram
    hfig_null = figure('Color','w','Position',[100 100 600 450]);
    hold on;
    histogram(xpos(maxIdx_perm), 44, ...
        'FaceColor',[0.75 0.75 0.75],'EdgeColor','k','LineWidth',0.8);
    xline(xpos(obsMaxIdx),'-','Color',[0.2 0 0],'LineWidth',2);
    xlabel('Cranio-caudal position (mm)','FontSize',14);
    ylabel('Count','FontSize',14);
    legend({'Null maxima','Observed maximum'},'Location','best','FontSize',14,'Box','off');
    set(gca,'FontSize',14,'LineWidth',1.2,'TickDir','out'); box off;
    title(sprintf('%s — null maxima (BS)', sub),'Interpreter','none');
    savefig(hfig_null, fullfile(fig_dir, ...
        sprintf('step2_sub%s_null_maxima%s.fig', sub, out_suffix)));
    close(hfig_null);

    %% -----------------------------------------------------------------------
    %  P1 ONLY: Null distribution stratified by dominant orientation axis
    %% -----------------------------------------------------------------------
    if strcmp(sub, 'OP00212')
        % For each permutation: find max source, get its filter, SVD ->
        % dominant orientation axis (x=1, y=2, z=3)
        perm_axis = zeros(nPerm,1);
        for ip = 1:nPerm
            src_idx = maxIdx_perm(ip);
            F = filters{src_idx};
            if isempty(F)
                perm_axis(ip) = 0;
                continue
            end
            [U,~,~]  = svd(F, 'econ');
            dom_vec  = U(:,1);
            [~, ax]  = max(abs(dom_vec));
            perm_axis(ip) = ax;
        end

        ax_colors = [0.22 0.49 0.82;   % x — right-left (blue)
                     0.86 0.37 0.13;   % y — cranial-caudal (orange)
                     0.18 0.63 0.27];  % z — dorsal-ventral (green)

        % Clamp negative maxPerm to 0 (one-sided test)
        maxPerm_pos = max(maxPerm, 0);
        nBins       = 30;
        bin_edges   = linspace(0, max(maxPerm_pos)*1.05, nBins+1);
        bin_ctrs    = (bin_edges(1:end-1) + bin_edges(2:end)) / 2;

        counts = zeros(nBins, 3);
        for ax_id = 1:3
            counts(:,ax_id) = histcounts(maxPerm_pos(perm_axis == ax_id), bin_edges)';
        end

        x_uniform  = (1:nBins)';
        tick_step  = max(1, floor(nBins/8));
        tick_idx   = 1:tick_step:nBins;
        scale_fac  = 1e4;

        hfig_nulldist = figure('Color','w','Position',[100 100 500 430]);
        hb = bar(x_uniform, counts, 'grouped', 'BarWidth', 0.85);
        for ax_id = 1:3
            hb(ax_id).FaceColor = ax_colors(ax_id,:);
            hb(ax_id).EdgeColor = 'none';
            hb(ax_id).FaceAlpha = 0.85;
        end
        set(gca, 'XTick', tick_idx, ...
            'XTickLabel', arrayfun(@(v) sprintf('%.2f', v*scale_fac), ...
            bin_ctrs(tick_idx), 'UniformOutput', false));
        xtickangle(35);
        hold on;
        [~, thr_bin] = min(abs(bin_ctrs - thr95));
        xline(thr_bin, '--k', 'LineWidth', 2, 'HandleVisibility','off');
        obs_max = max(coh_diff);
        [~, obs_bin] = min(abs(bin_ctrs - obs_max));
        xline(obs_bin, '-', 'Color',[0.55 0 0], 'LineWidth', 2.5, 'HandleVisibility','off');
        xlabel('Max coherence diff (static \minus rest)  \times 10^{-4}','FontSize',14);
        ylabel('Count','FontSize',14);
        set(gca,'FontSize',13,'LineWidth',1.1,'TickDir','out'); box off;
        legend off;
        title(sprintf('%s — null distribution by orientation (smoothed %d mm, BS)', ...
            sub, fwhm_mm),'Interpreter','none','FontSize',13);
        savefig(hfig_nulldist, fullfile(fig_dir, ...
            sprintf('step2_sub%s_null_distribution_by_orientation%s.fig', sub, out_suffix)));
        close(hfig_nulldist);
    end

    %% Interpolate onto spinal mesh
    source_p         = coh_source;
    source_p.avg.coh = invp_smooth;
    cfg_interp = []; cfg_interp.parameter = 'coh';
    spine_int = ft_sourceinterpolate(cfg_interp, source_p, mesh_wm);

    mesh_cut  = clip_torso(mesh_torso);
    invp_max  = max(invp_smooth);
    if invp_max <= invpthr
        clim_spine = [invpthr invpthr + 0.5];
    else
        clim_spine = [invpthr invp_max];
    end

    %% Figure — spine mesh with torso
    hfig_spine = figure;
    cfg_plot = []; cfg_plot.figure = 'gcf'; cfg_plot.method = 'surface';
    cfg_plot.funparameter = 'coh'; cfg_plot.funcolormap = cmap_hot;
    cfg_plot.funcolorlim = clim_spine; cfg_plot.projmethod = 'nearest';
    cfg_plot.surffile = mesh_wm;
    ft_sourceplot(cfg_plot, spine_int);
    view(-250,-1); camlight; ax = gca; ax.FontSize = 14;
    hold on;
    ft_plot_mesh(mesh_brain,'facecolor',[0.8 0.3 0.3],'facealpha',0.07,'edgecolor','none');
    ft_plot_mesh(mesh_cut,  'facecolor',[0.3 0.3 0.9],'facealpha',0.1, 'edgecolor','none');
    ft_plot_mesh(mesh_bone, 'facecolor',[0.9 0.85 0.7],'facealpha',0.3, 'edgecolor','none');
    ft_plot_sens(ftdat.grad,'coilshape','point','coilsize',6);
    ft_plot_mesh(mesh_lungs,'facecolor',[0.8 0.3 0.3],'facealpha',0.1, 'edgecolor','none');
    ft_plot_mesh(mesh_heart,'facecolor',[0.8 0.3 0.3],'facealpha',0.1, 'edgecolor','none');
    title(sprintf('%s — spine-EMG coherence (BS)', sub),'Interpreter','none');
    savefig(hfig_spine, fullfile(fig_dir, ...
        sprintf('step2_sub%s_spineEMG_coherence%s.fig', sub, out_suffix)));
    close(hfig_spine);

    %% Figure — mesh only for P1
    if ss == 1
        hfig_mesh = figure('Color','w');
        cfg_m = []; cfg_m.figure = 'gcf'; cfg_m.method = 'surface';
        cfg_m.funparameter = 'coh'; cfg_m.funcolormap = cmap_hot;
        cfg_m.funcolorlim = clim_spine; cfg_m.projmethod = 'nearest';
        cfg_m.surffile = mesh_wm;
        ft_sourceplot(cfg_m, spine_int);
        colorbar off;
        view(-250,-1); camlight; ax = gca; ax.FontSize = 14;
        title('Participant 1 — spine-EMG coherence (BS)', ...
            'Interpreter','none','FontSize',13);
        sig_invp = invp_smooth; sig_invp(~mask) = -inf;
        if any(isfinite(sig_invp))
            peaks_p1 = sources_cent.pos(sig_invp >= ...
                max(sig_invp(isfinite(sig_invp)))*0.99, :);
            hold on;
            scatter3(peaks_p1(:,1), peaks_p1(:,2), peaks_p1(:,3)+10, ...
                200, '.', 'filled', ...
                'MarkerFaceColor',[1 1 0],'MarkerEdgeColor','k','LineWidth',1.5);
            scatter_obj = findobj(gca,'Type','Scatter');
            uistack(scatter_obj,'top');
        end
        savefig(hfig_mesh, fullfile(fig_dir, ...
            sprintf('step2_sub%s_spineEMG_meshonly%s.fig', sub, out_suffix)));
        close(hfig_mesh);
    end

    %% Store results
    [~, obsIdx] = max(invp_smooth);
    subjResults(ss).subjID      = sub;
    subjResults(ss).coh_diff    = coh_diff;
    subjResults(ss).thr95       = thr95;
    subjResults(ss).sig_mask    = mask;
    subjResults(ss).pos         = sources_cent.pos;
    subjResults(ss).inside      = sources_cent.inside;
    subjResults(ss).maxdiff.idx = obsIdx;
    subjResults(ss).maxdiff.pos = sources_cent.pos(obsIdx,:);
    subjResults(ss).invp_smooth = invp_smooth;
    % dom_ori already set above

    % Save per-subject
    save(fullfile(save_dir, sprintf('subResult_sub%s%s.mat', sub, out_suffix)), ...
        'coh_diff','cohDiff_perm','thr95','mask','invp_smooth','pvals');
end

%% Save group results
save(fullfile(save_dir, ['groupRes_spine_DICS_' out_suffix '.mat']), 'subjResults');
fprintf('\nAll subjects complete. Running group analysis...\n');

%% =========================================================================
%  GROUP ANALYSIS
%% =========================================================================
run_group_spine(subjResults, sources_cent, save_dir, out_suffix, cmap_hot, ...
    mesh_wm, mesh_brain, mesh_torso, mesh_bone, mesh_lungs, mesh_heart, ...
    fwhm_mm, fig_dir);

fprintf('\n=== DONE ===\n');

%% =========================================================================
%  LOCAL FUNCTIONS
%% =========================================================================

function run_group_spine(subjResults, sources_cent, save_dir, out_suffix, cmap_hot, ...
                          mesh_wm, mesh_brain, mesh_torso, mesh_bone, mesh_lungs, mesh_heart, ...
                          fwhm, fig_dir)

nSubjects     = length(subjResults);
nsourcepoints = size(sources_cent.pos,1);
all_masks     = zeros(nsourcepoints, nSubjects);
for s = 1:nSubjects
    m = double(subjResults(s).coh_diff > subjResults(s).thr95);
    all_masks(:,s) = m(:);
end

sig_pos = false(nSubjects,1);
for s = 1:nSubjects
    sig_pos(s) = any(all_masks(:,s));
end
fprintf('  %d/%d subjects show significant spine-EMG coherence (smoothed %d mm, BS)\n', ...
    sum(sig_pos), nSubjects, fwhm);

threshold      = 0.2;
prevalence_loc = mean(all_masks, 2);

group_ft        = [];
group_ft.pos    = sources_cent.pos;
group_ft.inside = sources_cent.inside;
group_ft.pow    = prevalence_loc;
group_ft.pow(group_ft.pow < threshold) = 0;

cfg = []; cfg.parameter = 'pow'; cfg.interpmethod = 'nearest';
group_int = ft_sourceinterpolate(cfg, group_ft, mesh_wm);

mesh_cut = clip_torso(mesh_torso);

%% Group prevalence — full figure with torso
hfig_prev = figure;
cfg2 = []; cfg2.method = 'surface'; cfg2.funparameter = 'pow';
cfg2.maskparameter = 'mask';
cfg2.funcolorlim   = [threshold max(group_int.pow)];
cfg2.funcolormap   = cmap_hot;
cfg2.projmethod    = 'nearest';
cfg2.surffile      = mesh_wm;
cfg2.opacitylim    = [threshold max(group_int.pow)];
cfg2.opacitymap    = 'rampup';
ft_sourceplot(cfg2, group_int);
view(-250,-1); camlight; ax = gca; ax.FontSize = 14; hold on;
ft_plot_mesh(mesh_brain,'facecolor',[0.8 0.3 0.3],'facealpha',0.07,'edgecolor','none');
ft_plot_mesh(mesh_cut,  'facecolor',[0.3 0.3 0.9],'facealpha',0.1, 'edgecolor','none');
ft_plot_mesh(mesh_bone, 'facecolor',[0.9 0.85 0.7],'facealpha',0.3, 'edgecolor','none');
ft_plot_mesh(mesh_lungs,'facecolor',[0.8 0.3 0.3],'facealpha',0.1, 'edgecolor','none');
ft_plot_mesh(mesh_heart,'facecolor',[0.8 0.3 0.3],'facealpha',0.1, 'edgecolor','none');
title(sprintf('Group prevalence — spine-EMG (smoothed %d mm, BS)', fwhm), ...
    'Interpreter','none');
savefig(hfig_prev, fullfile(fig_dir, ['step2_group_spineEMG_prevalence' out_suffix '.fig']));
close(hfig_prev);

%% Group prevalence — mesh only with star at peak
hfig_grp_mesh = figure('Color','w','Position',[100 100 400 650]);
cfg3 = []; cfg3.method = 'surface'; cfg3.funparameter = 'pow';
cfg3.funcolorlim = [threshold max(group_int.pow)];
cfg3.funcolormap = cmap_hot;
cfg3.projmethod  = 'nearest';
cfg3.surffile    = mesh_wm;
cfg3.opacitylim  = [threshold max(group_int.pow)];
cfg3.opacitymap  = 'rampup';
ft_sourceplot(cfg3, group_int);
colorbar off;
view(-250,-1); camlight; ax = gca; ax.FontSize = 14;
title(sprintf('Group prevalence — spine-EMG (smoothed %d mm, BS)', fwhm), ...
    'Interpreter','none','FontSize',13);
prev_vals = group_ft.pow; prev_vals(prev_vals < threshold) = -inf;
if any(isfinite(prev_vals))
    peaks_grp = sources_cent.pos(prev_vals >= ...
        max(prev_vals(isfinite(prev_vals)))*0.99, :);
    hold on;
    scatter3(peaks_grp(:,1), peaks_grp(:,2), peaks_grp(:,3)+10, ...
        200, 'p', 'filled', ...
        'MarkerFaceColor',[1 1 0],'MarkerEdgeColor','k','LineWidth',1.5);
    scatter_obj = findobj(gca,'Type','Scatter');
    uistack(scatter_obj,'top');
end
savefig(hfig_grp_mesh, fullfile(fig_dir, ...
    ['step2_group_spineEMG_prevalence_meshonly' out_suffix '.fig']));
close(hfig_grp_mesh);

%% Subject line plot
subj_cmap = [27,158,119; 217,95,2; 117,112,179; 231,41,138; 102,166,30; ...
             230,171,2;  166,118,29; 102,102,102; 55,126,184] / 255;
x = sources_cent.pos(:,2);
hfig_lines = figure; hold on;
for s = 1:nSubjects
    cdiff = subjResults(s).coh_diff;
    thr   = subjResults(s).thr95;
    sig   = cdiff > thr;
    c     = subj_cmap(s,:);
    for i = 1:length(x)-1
        if sig(i) && sig(i+1)
            plot(x(i:i+1), cdiff(i:i+1), '-', 'Color', c, ...
                'LineWidth',1.5,'HandleVisibility','off');
        else
            plot(x(i:i+1), cdiff(i:i+1), '-', 'Color',[0.7 0.7 0.7], ...
                'LineWidth',1,'HandleVisibility','off');
        end
    end
    plot(x(sig), cdiff(sig), '.', 'Color', c, 'MarkerSize',12,'HandleVisibility','off');
    if sig_pos(s)
        h(s) = plot(nan, nan, '-', 'Color', c, 'LineWidth',1.5);
    else
        h(s) = plot(nan, nan, '-', 'Color',[0.7 0.7 0.7], 'LineWidth',1.5);
    end
end
yline(0,':k','HandleVisibility','off');
xlabel('Cranial-caudal position (mm)'); ylabel('Coherence difference (stat-rest)');
title(sprintf('Spine-EMG coherence differences (smoothed %d mm, BS)', fwhm));
legend(h, arrayfun(@(s) sprintf('Participant %d',s), 1:nSubjects,'UniformOutput',false), ...
    'Location','bestoutside');
set(gca,'FontSize',13); grid on;
savefig(hfig_lines, fullfile(fig_dir, ...
    ['step2_group_spineEMG_subject_lines' out_suffix '.fig']));
close(hfig_lines);

%% -----------------------------------------------------------------------
%  FIGURE C — Group dominant orientation bar chart
%  Colours match paper figure: purple=right-left(x), dark blue=cranial-
%  caudal(y), salmon=dorsal-ventral(z). P1 boxed.
%
%  Uses subjResults(s).dom_ori (computed once per subject in the main
%  loop, no need to recompute from coh_source/maxIdx_perm here, those
%  only exist in the subject-loop workspace and are not passed into this
%  function).
%% -----------------------------------------------------------------------
ax_colors = [0.55 0.20 0.75;   % x — right-left (purple)
             0.15 0.25 0.60;   % y — cranial-caudal (dark blue)
             0.90 0.40 0.30];  % z — dorsal-ventral (salmon)
ax_names  = {'Right-left','Cranial-caudal','Dorsal-ventral'};

ori_mat = nan(nSubjects, 3);
for s = 1:nSubjects
    if isfield(subjResults(s),'dom_ori') && ~isempty(subjResults(s).dom_ori) ...
            && ~any(isnan(subjResults(s).dom_ori))
        ori_mat(s,:) = subjResults(s).dom_ori(:)';
    end
end

hfig_ori = figure('Color','w','Position',[100 100 700 380]);
hold on;
bw      = 0.25;
offsets = [-bw, 0, bw];
for ax_id = 1:3
    for s = 1:nSubjects
        if isnan(ori_mat(s,ax_id)), continue; end
        b = bar(s + offsets(ax_id), ori_mat(s,ax_id), bw*0.9);
        b.FaceColor = ax_colors(ax_id,:);
        b.EdgeColor = 'none';
        b.FaceAlpha = 0.85;
    end
end
% Dummy bars for legend
for ax_id = 1:3
    bar(nan, nan, 'FaceColor', ax_colors(ax_id,:), 'EdgeColor','none', ...
        'DisplayName', ax_names{ax_id});
end
% Box around P1
rectangle('Position',[0.55 0 1 1.12], 'EdgeColor','k', 'LineWidth',2, ...
    'FaceColor','none');
set(gca, 'XTick', 1:nSubjects, ...
    'XTickLabel', arrayfun(@(s) sprintf('P%d',s), 1:nSubjects, 'UniformOutput',false));
xlabel('Participant','FontSize',13);
ylabel('Axis alignment (|component|, unit vector)','FontSize',13);
title('Spinal cord source orientation','FontSize',13,'Interpreter','none');
legend('Location','northeast','Box','off','FontSize',11);
ylim([0 1.15]);
set(gca,'FontSize',12); box off;
savefig(hfig_ori, fullfile(fig_dir, ['step2_group_dominant_orientation' out_suffix '.fig']));
saveas(hfig_ori,  fullfile(fig_dir, ['step2_group_dominant_orientation' out_suffix '.png']));
close(hfig_ori);

%% Cluster ROI for VE
mask_thresh = prevalence_loc >= threshold;
pos_thresh  = sources_cent.pos(mask_thresh,:);
if ~isempty(pos_thresh)
    distMat     = squareform(pdist(pos_thresh));
    G           = graph(distMat < 6);
    bins        = conncomp(G);
    [~, idxMax] = max(histcounts(bins, 1:(max(bins)+1)));
    ROIpos      = pos_thresh(bins == idxMax, :);

    hfig_clust = figure; hold on;
    plot(sources_cent.pos(:,2), prevalence_loc);
    for k = 1:size(ROIpos,1)
        plot(ROIpos(k,2), 0.2, 'r*');
    end
    xlabel('Cranial-caudal position (mm)'); ylabel('Prevalence');
    title('Spinal prevalence + VE cluster (BS)','Interpreter','none');
    savefig(hfig_clust, fullfile(fig_dir, ...
        ['step2_spinal_prevalence_VE_cluster' out_suffix '.fig']));
    close(hfig_clust);

    save(fullfile(save_dir, ['cluster_spineEMG_pos' out_suffix '.mat']), 'ROIpos');
    fprintf('  ROI cluster saved: %d sources\n', size(ROIpos,1));
else
    fprintf('  WARNING: no sources exceeded prevalence threshold — ROI cluster not saved\n');
end

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
function thr = compute_threshold(cohDiff_perm, mult_comp_corr, nsourcepoints)
    maxPerm = max(cohDiff_perm, [], 1);
    if mult_comp_corr
        thr = prctile(maxPerm, 95);
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
function invp_smooth = smooth_invp(coh_diff, cohDiff_perm, nsourcepoints, nPerm)
    invp_smooth = zeros(nsourcepoints,1);
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
function mesh_cut = clip_torso(mesh_torso)
    y = mesh_torso.vertices(:,2);
    keep_vert = y > -200;
    new_idx   = zeros(size(keep_vert));
    new_idx(keep_vert) = 1:sum(keep_vert);
    faces_keep        = all(keep_vert(mesh_torso.faces), 2);
    mesh_cut.vertices = mesh_torso.vertices(keep_vert,:);
    mesh_cut.faces    = new_idx(mesh_torso.faces(faces_keep,:));
    mesh_cut.unit     = mesh_torso.unit;
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
