%% step2_spineEMG_coherence_DICS.m
% Step 2: Spine-EMG DICS coherence with Biot-Savart (BS) forward model.
%
% Uses the BS leadfield with label-based channel matching. Lambda = 10%,
% smoothing 20mm FWHM. All 9 spine subjects. Outputs saved with _BS
% suffix to distinguish from older (BEM) results.
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
%   subResult_sub<ID>_BS.mat                          - per-subject coh_diff, thr95, mask, invp_smooth, pvals
%   groupRes_spine_DICS__BS.mat                        - subjResults struct for all subjects
%   cluster_spineEMG_pos_BS.mat                        - VE ROI cluster positions
%   step2_sub<ID>_spineEMG_coherence_BS.fig            - per-subject coherence map (with torso)
%   step2_sub1_spineEMG_meshonly_BS.fig                - Participant 1 mesh-only map
%   step2_sub<ID>_orientation_BS.fig                   - per-subject source orientation (null vs observed)
%   step2_group_spineEMG_prevalence_BS.fig             - group prevalence map (with torso)
%   step2_group_spineEMG_prevalence_meshonly_BS.fig    - group prevalence, mesh only
%   step2_group_spineEMG_subject_lines_BS.fig          - per-subject coherence-difference traces
%   step2_spinal_prevalence_VE_cluster_BS.fig          - prevalence + VE cluster ROI
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
lf_path = 'C:\Leadfields meshes\leadfield_experimental_bslaw_experimental.mat';

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
%  - geometries_experimental: all spine/torso meshes + sources_cent
%  - geometries_cervical_realistic: mesh_brain only (other meshes outdated)
%% =========================================================================
fprintf('Loading geometry...\n');
geom_exp = load(geomfile);

sources_cent = geom_exp.sources_cent;
mesh_torso   = geom_exp.mesh_torso;
mesh_wm      = geom_exp.mesh_wm;
mesh_bone    = geom_exp.mesh_bone;
mesh_lungs   = geom_exp.mesh_lungs;
mesh_heart   = geom_exp.mesh_heart;

% Load mesh_brain from old geomfile — only this mesh, nothing else
warning('Need to add brain mesh to new geom!!!!')
geom_brain = load(geomfile_brain, 'mesh_brain');
mesh_brain = geom_brain.mesh_brain;

mesh_wm.unit = 'mm';

nsourcepoints = size(sources_cent.pos, 1);
fprintf('  Source space: %d points\n', nsourcepoints);

%% =========================================================================
%  BUILD SMOOTHER ONCE (shared across all subjects)
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
    % keeptrials='yes' for stat/rest (permutation test needs trial-level)
    % keeptrials='no' for combined (common spatial filter)
    cfg_fr = []; cfg_fr.output = 'powandcsd'; cfg_fr.method = 'mtmfft';
    cfg_fr.foilim = fband; cfg_fr.tapsmofrq = 1; cfg_fr.keeptrials = 'yes';
    cfg_av = []; cfg_av.avgoverfreq = 'yes';
    cfg_sel = []; cfg_sel.channel = [Lf.label; {'EXG1'}];

    freqdat_tr = ft_freqanalysis(cfg_fr, ftdat);
    freqdat_tr = ft_selectdata(cfg_av,  freqdat_tr);

    trialinfo = ftdat.trialinfo;
    statidx   = find(trialinfo == 1);
    restidx   = find(trialinfo == 2);
    nTrials   = min(numel(statidx), numel(restidx));
    fprintf('  nTrials (stat/rest): %d\n', nTrials);

    cfg = []; cfg.trials = statidx(1:nTrials);
    statdat = ft_selectdata(cfg, freqdat_tr);
    cfg = []; cfg.trials = restidx(1:nTrials);
    restdat = ft_selectdata(cfg, freqdat_tr);

    % Combined for common spatial filter
    cfg_fr2 = []; cfg_fr2.output = 'powandcsd'; cfg_fr2.method = 'mtmfft';
    cfg_fr2.foilim = fband; cfg_fr2.tapsmofrq = 1; cfg_fr2.keeptrials = 'no';
    freqdat = ft_freqanalysis(cfg_fr2, ftdat);
    cfg_av2 = []; cfg_av2.avgoverfreq = 'yes';
    freqdat = ft_selectdata(cfg_av2, freqdat);

    % Select channels
    statdat = ft_selectdata(cfg_sel, statdat);
    restdat = ft_selectdata(cfg_sel, restdat);
    freqdat = ft_selectdata(cfg_sel, freqdat);

    %% Sourcemodel
    sourcemodel = [];
    sourcemodel.pos       = Lf.pos;
    sourcemodel.unit      = 'mm';
    sourcemodel.inside    = logical(Lf.inside);
    sourcemodel.leadfield = Lf.leadfield;
    sourcemodel.label     = Lf.label;

    %% Common spatial filter from combined data
    fprintf('  Computing common spatial filter...\n');
    cfg_dics = [];
    cfg_dics.sourcemodel     = sourcemodel;
    cfg_dics.headmodel       = dummyvol;
    cfg_dics.dics.keepfilter = 'yes';
    cfg_dics.dics.lambda     = lambda;
    cfg_dics.method          = 'dics';
    cfg_dics.refchan         = 'EXG1';
    coh_source = ft_sourceanalysis(cfg_dics, freqdat);

    %% Permutation test — contraction vs rest
    fprintf('  Running permutation test (%d permutations)...\n', numpermutation);
    cfg_perm = [];
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

    %% Null maxima orientation diagnostic — all three components
    [~, maxIdx_perm] = max(cohDiff_perm, [], 1);

    filters     = coh_source.avg.filter;
    axis_labels = {'RL (x)', 'CC (y)', 'VD (z)'};
    ax_colors   = [0.22 0.49 0.82;
                   0.86 0.37 0.13;
                   0.18 0.63 0.27];

    % For each permutation: extract all 3 SVD component magnitudes at peak location
    perm_comps = nan(nPerm, 3);
    for ip = 1:nPerm
        src_idx = maxIdx_perm(ip);
        F = filters{src_idx};
        if isempty(F), continue; end
        [U,~,~] = svd(F, 'econ');
        dom_vec = U(:,1);
        perm_comps(ip,:) = abs(dom_vec) / sum(abs(dom_vec));  % normalise to sum=1
    end

    % Observed orientation at peak coherence location
    [~, obsMaxIdx] = max(coh_diff);
    F_obs = filters{obsMaxIdx};
    [U_obs,~,~] = svd(F_obs, 'econ');
    dom_obs      = U_obs(:,1);
    obs_comps    = abs(dom_obs) / sum(abs(dom_obs));

    % Figure — dot plot: observed vs null median ± IQR
    hfig_ori = figure('Color','w','Position',[100 100 420 400]);
    hold on;

    for ax_id = 1:3
        vals     = perm_comps(:, ax_id);
        vals     = vals(~isnan(vals));
        null_med = median(vals);
        null_iqr = iqr(vals);
        col      = ax_colors(ax_id,:);

        % Null median ± IQR/2
        errorbar(ax_id, null_med, null_iqr/2, 's', ...
            'Color', col, 'MarkerFaceColor', col, ...
            'MarkerSize', 10, 'LineWidth', 2.0, ...
            'CapSize', 8, 'HandleVisibility','off');

        % Observed value
        scatter(ax_id + 0.18, obs_comps(ax_id), 60, [0.4 0.4 0.4], 'filled', ...
            'MarkerEdgeColor','none', 'MarkerFaceAlpha', 0.5, ...
            'HandleVisibility','off');
    end

    % Legend
    h_null = plot(nan, nan, 's', 'Color',[0.5 0.5 0.5], ...
        'MarkerFaceColor',[0.5 0.5 0.5], 'MarkerSize',10, 'LineWidth',2);
    h_obs = scatter(nan, nan, 60, [0.4 0.4 0.4], 'filled', ...
        'MarkerEdgeColor','none', 'MarkerFaceAlpha', 0.5);
    legend([h_null h_obs], {'Null median ± IQR','Observed'}, ...
        'Location','northwest','Box','off','FontSize',12);

    set(gca,'XTick',1:3,'XTickLabel',axis_labels,'FontSize',13);
    ylabel('Normalised component magnitude','FontSize',13);
    ylim([0 1]);
    title(sprintf('%s — source orientation', sub), ...
        'Interpreter','none','FontSize',12,'FontWeight','normal');
    box off;
    savefig(hfig_ori, fullfile(fig_dir, ...
        sprintf('step2_sub%s_orientation%s.fig', sub, out_suffix)));
    close(hfig_ori);

    %% Interpolate onto spinal mesh
    source_p = coh_source;
    source_p.avg.coh = invp_smooth;
    cfg_interp = []; cfg_interp.parameter = 'coh';
    spine_int = ft_sourceinterpolate(cfg_interp, source_p, mesh_wm);

    mesh_cut = clip_torso(mesh_torso);

    invp_max = max(invp_smooth);
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

    %% Figure — mesh only for subject 1
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
        % Star at peak significant source
        sig_invp = invp_smooth; sig_invp(~mask) = -inf;
        if any(isfinite(sig_invp))
            peaks_p1 = sources_cent.pos(sig_invp >= ...
                max(sig_invp(isfinite(sig_invp)))*0.99, :);
            hold on;
            scatter3(peaks_p1(:,1), peaks_p1(:,2), peaks_p1(:,3)+10, ...
                200, '.', 'filled', ...
                'MarkerFaceColor',[1 1 0], 'MarkerEdgeColor','k','LineWidth',1.5);
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

    % Save per-subject immediately in case of crash
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

group_ft         = [];
group_ft.pos     = sources_cent.pos;
group_ft.inside  = sources_cent.inside;
group_ft.pow     = prevalence_loc;
group_ft.pow(group_ft.pow < threshold) = 0;

cfg = []; cfg.parameter = 'pow'; cfg.interpmethod = 'nearest';
group_int = ft_sourceinterpolate(cfg, group_ft, mesh_wm);

mesh_cut = clip_torso(mesh_torso);

% Group prevalence — full figure with torso
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

% Group prevalence — mesh only with star at peak
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
        'MarkerFaceColor',[1 1 0], 'MarkerEdgeColor','k','LineWidth',1.5);
    scatter_obj = findobj(gca,'Type','Scatter');
    uistack(scatter_obj,'top');
end
savefig(hfig_grp_mesh, fullfile(fig_dir, ...
    ['step2_group_spineEMG_prevalence_meshonly' out_suffix '.fig']));
close(hfig_grp_mesh);

% Subject line plot
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

% Cluster ROI for VE
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
