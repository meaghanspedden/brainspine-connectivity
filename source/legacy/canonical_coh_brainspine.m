%% step_cancoh_manual.m
% Manual implementation of canonical coherence (caCOH) between M1 VE
% and spinal cord source points, following Vidaurre et al.
% Since M1 is a single channel, w_B is scalar and cancels — we only
% need to find the cord weight vector w_A that maximises coherence with M1.
%
% Algorithm:
%   1. Compute CSD matrix: S_AA (53x53 cord), S_AB (53x1 cord-M1), S_BB (scalar)
%   2. SVD dimensionality reduction of Re(S_AA) — keep 99% variance
%   3. Search over phase phi in [0,pi): rotate S_AB, solve eigenvalue problem
%   4. At optimal phi: extract w_A and spatial pattern a_A = Re(S_AA) * w_A
%   5. Map a_A onto cord mesh and compare with spine-EMG prevalence map

clear all; close all; clc;

%% =========================================================================
%  USER CONFIG
%% =========================================================================
fieldtrip_path = 'C:\Users\mspedden\Documents\fieldtrip';
spm_path       = 'C:\Users\mspedden\Documents\spm';
bsc_path       = 'C:\Users\mspedden\Documents\brainspineconnectivity\source';
data_root      = 'C:\spinecoh_data';
save_dir       = 'C:\Users\mspedden\Documents\brainspine_savetest';
geomfile_spine = 'C:\Leadfields meshes\geometries_experimental.mat';
geomfile_brain = 'C:\Leadfields meshes\geometries_cervical_realistic.mat';
lf_path        = 'C:\Leadfields meshes\leadfield_experimental_bslaw_experimental.mat';

subs_brain = {'OP00212','OP00213','OP00215','OP00219', ...
              'OP00225','OP00221','OP00224'};
subs_spine = {'OP00212','OP00213','OP00215','OP00219', ...
              'OP00220','OP00221','OP00224','OP00225','OP00226'};
subs_both  = intersect(subs_brain, subs_spine, 'stable');

brain_suffix  = '_brain_pct10';
spine_suffix  = '_bslaw';
out_suffix    = '_cancoh_manual_bslaw';

fband         = [10 35];
n_phi         = 50;        % phase search resolution
svd_var_thr   = 0.90;      % keep components explaining 99% variance
reg_alpha     = 0.1;      % Tikhonov regularisation (5% of max eigenvalue)

%% =========================================================================
%  SETUP
%% =========================================================================
addpath(bsc_path); addpath(spm_path);
spm('defaults','EEG');
addpath(fieldtrip_path);
ft_defaults;

if ~exist(save_dir,'dir'), mkdir(save_dir); end
fig_dir = fullfile(save_dir, 'figures', 'cancoh_manual');
if ~exist(fig_dir,'dir'), mkdir(fig_dir); end

ncol     = 256;
addpath(fullfile(fieldtrip_path,'external','matplotlib'));
cmap_hot = [[0.92 0.92 0.92]; flipud(magma(ncol-1))];

%% =========================================================================
%  LOAD GEOMETRY
%% =========================================================================
fprintf('Loading geometry...\n');
geom_exp      = load(geomfile_spine);
sources_cent  = geom_exp.sources_cent;
mesh_torso    = geom_exp.mesh_torso;
mesh_wm       = geom_exp.mesh_wm;
mesh_wm.unit  = 'mm';

geom_brain    = load(geomfile_brain, 'mesh_brain');
mesh_brain    = geom_brain.mesh_brain;

nsourcepoints = size(sources_cent.pos, 1);
cord_pos      = sources_cent.pos(:,2);
fprintf('  Cord source points: %d\n', nsourcepoints);

%% =========================================================================
%  LOAD SPINE-EMG PREVALENCE MAP FOR OVERLAP COMPARISON
%% =========================================================================
spineEMG_file = fullfile(save_dir, ...
    ['groupRes_spine_DICS_bemv2' spine_suffix '.mat']);
if exist(spineEMG_file,'file')
    loaded_spine = load(spineEMG_file, 'subjResults');
    spine_subj   = loaded_spine.subjResults;
    nSpine       = numel(spine_subj);
    all_spine_masks = zeros(nsourcepoints, nSpine);
    for s = 1:nSpine
        all_spine_masks(:,s) = double(spine_subj(s).coh_diff > spine_subj(s).thr95);
    end
    spine_prevalence = mean(all_spine_masks, 2);
    fprintf('  Spine-EMG prevalence map loaded (%d subjects)\n', nSpine);
else
    spine_prevalence = [];
    fprintf('  WARNING: spine-EMG prevalence map not found — overlap plot skipped\n');
end

%% =========================================================================
%  SUBJECT LOOP
%% =========================================================================
subjResults = struct();

for ss = 1:numel(subs_both)
    sub = subs_both{ss};
    fprintf('\n========================================\n');
    fprintf('  Subject %s (%d/%d)\n', sub, ss, numel(subs_both));
    fprintf('========================================\n');
    t_sub = tic;

    %% Load M1 VE
    ve_file = fullfile(save_dir, ...
        sprintf('sub%s_VE_brain_M1%s.mat', sub, brain_suffix));
    assert(exist(ve_file,'file')==2, 'Brain VE not found: %s', ve_file);
    loaded_ve = load(ve_file, 'VE_brain');
    VE_brain  = loaded_ve.VE_brain;
    VE_brain.label{1} = 'M1';
    fprintf('  M1 VE loaded\n');

    %% Load spine data
    run = '001'; if strcmp(sub,'OP00224'), run = '002'; end
    datafile = fullfile(data_root, ['sub-' sub], 'ses-001', 'meg', ...
        sprintf('pmergedoe1000mspddfflo45hi45hfcstatic_%s_array1.mat', run));

    D       = spm_eeg_load(datafile);
    grad_mm = D.sensors('MEG');
    ftdat   = spm2fieldtrip(D);

    badchans = D.chanlabels(D.badchannels);
    cfg = []; cfg.channel = setdiff(ftdat.label, badchans);
    ftdat = ft_selectdata(cfg, ftdat);

    % Contraction trials only
    statidx = find(ftdat.trialinfo == 1);
    cfg = []; cfg.trials = statidx;
    ftdat_stat   = ft_selectdata(cfg, ftdat);
    cfg.trials   = statidx;
    VE_stat      = ft_selectdata(cfg, VE_brain);

    samp_rate = ftdat.hdr.Fs;

    %% Build BSLaw spine leadfield
    fprintf('  Building BSLaw spine leadfield...\n');
    lf_data = load(lf_path);
    lf_raw  = lf_data.leadfield_bs;

    data_meg_labels        = ftdat.label(~strcmp(ftdat.label,'EXG1'));
    [common_labels,idx_lf] = intersect(lf_raw.label, data_meg_labels, 'stable');
    fprintf('  Data MEG: %d  |  LF: %d  |  Matched: %d\n', ...
        numel(data_meg_labels), numel(lf_raw.label), numel(common_labels));

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

    %% LCMV beamformer — all cord sources
    fprintf('  Computing LCMV beamformer...\n');
    cfg_tl = []; cfg_tl.covariance = 'yes';
    cfg_tl.covariancewindow = 'all';
    cfg_tl.channel = Lf.label;
    tlock = ft_timelockanalysis(cfg_tl, ftdat);

    cfg_lcmv = [];
    cfg_lcmv.method                   = 'lcmv';
    cfg_lcmv.headmodel                = dummyvol;
    cfg_lcmv.sourcemodel.pos          = Lf.pos;
    cfg_lcmv.sourcemodel.unit         = 'mm';
    cfg_lcmv.sourcemodel.inside       = logical(Lf.inside);
    cfg_lcmv.sourcemodel.leadfield    = Lf.leadfield;
    cfg_lcmv.sourcemodel.label        = Lf.label;
    cfg_lcmv.lcmv.keepfilter          = 'yes';
    source_lcmv = ft_sourceanalysis(cfg_lcmv, tlock);

    %% Extract VE timeseries for all cord source points
    fprintf('  Extracting cord source timeseries...\n');
    cord_ts_all = [];   % nsources x nSamples

    cfg_sel = []; cfg_sel.channel = Lf.label;
    ftdat_meg = ft_selectdata(cfg_sel, ftdat_stat);

    for si = 1:nsourcepoints
        cfg_ve = [];
        cfg_ve.pos          = sources_cent.pos(si,:);
        cfg_ve.radius       = 1;
        cfg_ve.method       = 'svd';
        cfg_ve.numcomponent = 1;
        VE_cord = ft_virtualchannel(cfg_ve, ftdat_meg, source_lcmv);
        cord_ts_all = [cord_ts_all; [VE_cord.trial{:}]];
    end

    % M1 timeseries
    M1_ts = [VE_stat.trial{:}];   % 1 x nSamples

    fprintf('  Cord timeseries: %d x %d\n', size(cord_ts_all));

    %% Compute CSD over beta band using Welch method
    fprintf('  Computing CSD matrix over beta band...\n');
    [S_AA, S_AB, S_BB, freq_vec] = compute_csd_beta(cord_ts_all, M1_ts, ...
        samp_rate, fband);

    fprintf('  CSD computed: S_AA=%dx%d  S_AB=%dx1  S_BB=scalar\n', ...
        size(S_AA,1), size(S_AA,2), size(S_AB,1));

    %% SVD dimensionality reduction of Re(S_AA)
    fprintf('  SVD dimensionality reduction...\n');
    ReS_AA = real(S_AA);
    [U_svd, D_svd, ~] = svd(ReS_AA);
    sv = diag(D_svd);
    cum_var = cumsum(sv) / sum(sv);
    n_comp  = find(cum_var >= svd_var_thr, 1, 'first');
    fprintf('  Keeping %d/%d components (%.1f%% variance)\n', ...
        n_comp, nsourcepoints, 100*cum_var(n_comp));

    % Projection matrix
    P = U_svd(:, 1:n_comp);   % nsources x n_comp

    fprintf('  Singular values: '); fprintf('%.3f ', sv(1:min(10,end))'); fprintf('\n');
fprintf('  Cumulative variance at each component:\n');
for k = [1 2 3 4 5 10 20 30 nsourcepoints]
    if k <= nsourcepoints
        fprintf('    %d components: %.1f%%\n', k, 100*cum_var(k));
    end
end

    % Project CSD into reduced space
    S_AA_r = P' * S_AA * P;   % n_comp x n_comp
    S_AB_r = P' * S_AB;       % n_comp x 1

    %% Canonical coherence: search over phase phi
    fprintf('  Phase search for optimal canonical coherence...\n');
    phi_grid = linspace(0, pi, n_phi+1);
    phi_grid = phi_grid(1:end-1);   % exclude pi (same as 0)

    coh_vals = nan(n_phi, 1);
    wA_vals  = nan(n_comp, n_phi);

    ReS_AA_r = real(S_AA_r);

    % Regularise: Tikhonov with reg_alpha * max_eigenvalue
    ev_max   = max(eig(ReS_AA_r));
    reg_mat  = reg_alpha * ev_max * eye(n_comp);
    S_AA_reg = ReS_AA_r + reg_mat;

    for pi_idx = 1:n_phi
        phi = phi_grid(pi_idx);

        % Rotate cross-spectrum
        S_AB_rot = real(exp(-1i*phi) * S_AB_r);   % n_comp x 1

        % Solve for w_A: eigenvalue problem reduces to linear solve
        % since M1 is 1D, w_B is scalar => w_A = S_AA_reg \ S_AB_rot
        w_A = S_AA_reg \ S_AB_rot;

        % Compute canonical coherence
        num   = abs(w_A' * S_AB_r)^2;
        denom = real(w_A' * S_AA_r * w_A) * real(S_BB);
        if denom > 0
            coh_vals(pi_idx) = num / denom;
        end
        wA_vals(:, pi_idx) = w_A;
    end

    %% Refine with Levenberg-Marquardt style iterations (eq. 10 iterations)
    [~, best_idx] = max(coh_vals);
    phi_opt = phi_grid(best_idx);
    w_A_opt = wA_vals(:, best_idx);

    for iter = 1:10
        % Numerical first and second derivative
        dphi = 1e-4;
        coh_p = coh_from_phi(phi_opt+dphi, S_AA_reg, S_AA_r, S_AB_r, S_BB);
        coh_m = coh_from_phi(phi_opt-dphi, S_AA_reg, S_AA_r, S_AB_r, S_BB);
        coh_c = coh_from_phi(phi_opt,      S_AA_reg, S_AA_r, S_AB_r, S_BB);

        d1 = (coh_p - coh_m) / (2*dphi);
        d2 = (coh_p - 2*coh_c + coh_m) / dphi^2;

        if d2 < 0
            step = -d1/d2;
            step = max(-0.1, min(0.1, step));   % clip step size
            phi_opt = phi_opt + step;
            phi_opt = mod(phi_opt, pi);
        else
            break;
        end
    end

    % Final w_A at optimal phi
    S_AB_rot_opt = real(exp(-1i*phi_opt) * S_AB_r);
    w_A_opt      = S_AA_reg \ S_AB_rot_opt;

    % Final canonical coherence value
    coh_final = coh_from_phi(phi_opt, S_AA_reg, S_AA_r, S_AB_r, S_BB);
    fprintf('  Optimal phi: %.3f rad  |  Canonical coherence: %.4f\n', ...
        phi_opt, coh_final);

    %% Project weights back to full source space
    w_A_full = P * w_A_opt;   % nsources x 1

    %% Spatial pattern (topography) — eq. 14: a_A = Re(S_AA) * w_A
    a_A = real(S_AA) * w_A_full;   % nsources x 1

    % Normalise to [0 1] for visualisation
    a_A_norm = (a_A - min(a_A)) / (max(a_A) - min(a_A));

    fprintf('  Peak weight at source %d (y=%.1f mm)\n', ...
        find(a_A_norm==max(a_A_norm)), ...
        sources_cent.pos(a_A_norm==max(a_A_norm),2));

    %% Interpolate spatial pattern onto spine mesh
    source_p         = [];
    source_p.pos     = sources_cent.pos;
    source_p.inside  = true(nsourcepoints,1);
    source_p.avg.coh = a_A_norm;

    cfg_interp = []; cfg_interp.parameter = 'coh';
    spine_int  = ft_sourceinterpolate(cfg_interp, source_p, mesh_wm);

    mesh_cut = clip_torso(mesh_torso);
    clim = [0 1];

    %% Figure — spatial pattern on spine mesh
    hfig_spine = figure;
    cfg_plot = []; cfg_plot.figure = 'gcf'; cfg_plot.method = 'surface';
    cfg_plot.funparameter = 'coh'; cfg_plot.funcolormap = cmap_hot;
    cfg_plot.funcolorlim  = clim; cfg_plot.projmethod = 'nearest';
    cfg_plot.surffile = mesh_wm;
    ft_sourceplot(cfg_plot, spine_int);
    view(-250,-1); camlight; ax = gca; ax.FontSize = 14;
    hold on;
    ft_plot_mesh(mesh_brain,'facecolor',[0.8 0.3 0.3],'facealpha',0.07,'edgecolor','none');
    ft_plot_mesh(mesh_cut,  'facecolor',[0.3 0.3 0.9],'facealpha',0.1, 'edgecolor','none');
    title(sprintf('%s — caCOH spatial pattern (BSLaw)', sub),'Interpreter','none');
    savefig(hfig_spine, fullfile(fig_dir, ...
        sprintf('cancoh_sub%s_pattern%s.fig', sub, out_suffix)));
    close(hfig_spine);

    %% Store
    subjResults(ss).subjID    = sub;
    subjResults(ss).a_A_norm  = a_A_norm;
    subjResults(ss).a_A       = a_A;
    subjResults(ss).w_A_full  = w_A_full;
    subjResults(ss).coh_final = coh_final;
    subjResults(ss).phi_opt   = phi_opt;
    subjResults(ss).n_comp    = n_comp;

    save(fullfile(save_dir, sprintf('cancoh_result_sub%s%s.mat', sub, out_suffix)), ...
        'a_A_norm','a_A','w_A_full','coh_final','phi_opt','n_comp');

    fprintf('  Subject time: %.1f min\n', toc(t_sub)/60);
end

%% Save group
save(fullfile(save_dir, ['groupRes_cancoh' out_suffix '.mat']), 'subjResults');
fprintf('\nAll subjects complete. Running group analysis...\n');

%% =========================================================================
%  GROUP ANALYSIS
%% =========================================================================
nSubjects = numel(subjResults);

%% Group mean spatial pattern
all_patterns = cat(2, subjResults(:).a_A_norm);
mean_pattern = mean(all_patterns, 2);
mean_pattern_norm = (mean_pattern - min(mean_pattern)) / ...
                    (max(mean_pattern) - min(mean_pattern));

%% Overlap with spine-EMG prevalence
if ~isempty(spine_prevalence)
    [r_overlap, p_overlap] = corr(mean_pattern_norm, spine_prevalence, ...
        'type','Spearman');
    fprintf('  Overlap with spine-EMG prevalence: r=%.3f p=%.3f\n', ...
        r_overlap, p_overlap);
end

%% Group mean pattern on mesh
source_grp         = [];
source_grp.pos     = sources_cent.pos;
source_grp.inside  = true(nsourcepoints,1);
source_grp.avg.coh = mean_pattern_norm;

cfg_interp = []; cfg_interp.parameter = 'coh';
grp_int = ft_sourceinterpolate(cfg_interp, source_grp, mesh_wm);

hfig_grp = figure;
cfg_plot = []; cfg_plot.figure = 'gcf'; cfg_plot.method = 'surface';
cfg_plot.funparameter = 'coh'; cfg_plot.funcolormap = cmap_hot;
cfg_plot.funcolorlim  = [0 1]; cfg_plot.projmethod = 'nearest';
cfg_plot.surffile = mesh_wm;
ft_sourceplot(cfg_plot, grp_int);
view(-250,-1); camlight; ax = gca; ax.FontSize = 14;
hold on;
ft_plot_mesh(mesh_brain,'facecolor',[0.8 0.3 0.3],'facealpha',0.07,'edgecolor','none');
ft_plot_mesh(mesh_cut,  'facecolor',[0.3 0.3 0.9],'facealpha',0.1, 'edgecolor','none');
title('Group mean — caCOH spatial pattern (BSLaw)','Interpreter','none');
savefig(hfig_grp, fullfile(fig_dir, ['cancoh_group_pattern' out_suffix '.fig']));

%% Overlap figure — caCOH pattern vs spine-EMG prevalence
if ~isempty(spine_prevalence)
    hfig_overlap = figure('Color','w','Position',[100 100 600 500]);
    hold on;

    scatter(spine_prevalence, mean_pattern_norm, 80, ...
        [0.2 0.4 0.8], 'filled', 'MarkerFaceAlpha',0.7);

    % Highlight anatomical ROI points
    roi_idx = 14:27;
    scatter(spine_prevalence(roi_idx), mean_pattern_norm(roi_idx), 120, ...
        [0.8 0.2 0.2], 'filled', 'MarkerFaceAlpha',0.9, ...
        'DisplayName','C6-T1 ROI');

    % Regression line
    ok = isfinite(spine_prevalence) & isfinite(mean_pattern_norm);
    p_fit = polyfit(spine_prevalence(ok), mean_pattern_norm(ok), 1);
    x_fit = linspace(min(spine_prevalence), max(spine_prevalence), 100);
    plot(x_fit, polyval(p_fit,x_fit), 'k--', 'LineWidth',1.5);

    text(0.05, 0.92, sprintf('r_s=%.2f, p=%.3f', r_overlap, p_overlap), ...
        'Units','normalized','FontSize',12);

    xlabel('Spine-EMG coherence prevalence','FontSize',13);
    ylabel('caCOH spatial pattern (normalised)','FontSize',13);
    title('caCOH pattern vs Spine-EMG prevalence','Interpreter','none');
    legend('All sources','C6-T1 ROI','Location','best','Box','off','FontSize',11);
    set(gca,'FontSize',13); box on; grid on;
    savefig(hfig_overlap, fullfile(fig_dir, ...
        ['cancoh_vs_spineEMG_overlap' out_suffix '.fig']));
    saveas(hfig_overlap, fullfile(fig_dir, ...
        ['cancoh_vs_spineEMG_overlap' out_suffix '.png']));
end

%% Line plot along cord
subj_cmap = [27,158,119; 217,95,2; 117,112,179; 231,41,138; 102,166,30; ...
             230,171,2;  166,118,29] / 255;

hfig_line = figure('Color','w','Position',[100 100 700 400]);
hold on;
for s = 1:nSubjects
    if s == 1
        plot(cord_pos, subjResults(s).a_A_norm, 'k-', 'LineWidth',2.5);
    else
        plot(cord_pos, subjResults(s).a_A_norm, '-', ...
            'Color',[subj_cmap(s,:) 0.5], 'LineWidth',1.2);
    end
end
plot(cord_pos, mean_pattern_norm, 'k--', 'LineWidth',2, ...
    'DisplayName','Group mean');
if ~isempty(spine_prevalence)
    yyaxis right;
    plot(cord_pos, spine_prevalence, 'b:', 'LineWidth',1.5, ...
        'DisplayName','Spine-EMG prevalence');
    ylabel('Spine-EMG prevalence','FontSize',13);
    yyaxis left;
end
xlabel('Cranio-caudal position (mm)','FontSize',13);
ylabel('caCOH pattern weight (normalised)','FontSize',13);
title('caCOH spatial pattern along cord vs spine-EMG prevalence (BSLaw)', ...
    'Interpreter','none');
legend('Location','best','FontSize',10,'Box','off');
set(gca,'FontSize',13); grid on; box on;
savefig(hfig_line, fullfile(fig_dir, ['cancoh_group_cord_profile' out_suffix '.fig']));
saveas(hfig_line,  fullfile(fig_dir, ['cancoh_group_cord_profile' out_suffix '.png']));

fprintf('\n=== CANCOH DONE ===\n');
fprintf('  Canonical coherence values per subject:\n');
for s = 1:nSubjects
    fprintf('    %s: coh=%.4f  n_comp=%d\n', ...
        subjResults(s).subjID, subjResults(s).coh_final, subjResults(s).n_comp);
end

%% =========================================================================
%  LOCAL FUNCTIONS
%% =========================================================================
function [S_AA, S_AB, S_BB, freq_vec] = compute_csd_beta(cord_ts, M1_ts, ...
        samp_rate, fband)
    % Compute averaged CSD over beta band using Welch method
    seg_len  = 2 * samp_rate;   % 2 second segments
    n_fft    = seg_len;
    n_overlap = round(seg_len * 0.5);
    win      = hann(seg_len);
    nSamp    = size(cord_ts, 2);
    nSources = size(cord_ts, 1);
    nFreqs   = n_fft/2 + 1;

    freq_vec = (0:nFreqs-1) * samp_rate / n_fft;
    beta_mask = freq_vec >= fband(1) & freq_vec <= fband(2);

    % Accumulate CSD
    S_AA_all = zeros(nSources, nSources, sum(beta_mask));
    S_AB_all = zeros(nSources, 1,        sum(beta_mask));
    S_BB_all = zeros(1,        1,        sum(beta_mask));

    n_segs = 0;
    start_idx = 1;
    while start_idx + seg_len - 1 <= nSamp
        idx = start_idx : start_idx + seg_len - 1;

        % Apply window
        cord_seg = cord_ts(:,idx) .* win';
        M1_seg   = M1_ts(idx)    .* win';

        % FFT
        cord_fft = fft(cord_seg, n_fft, 2);
        M1_fft   = fft(M1_seg,   n_fft);

        % One-sided
        cord_fft = cord_fft(:, 1:nFreqs);
        M1_fft   = M1_fft(1:nFreqs);

        % Accumulate beta band
        cord_beta = cord_fft(:, beta_mask);   % nSources x nBetaFreqs
        M1_beta   = M1_fft(beta_mask);        % 1 x nBetaFreqs

        for fi = 1:sum(beta_mask)
            c  = cord_beta(:,fi);
            m  = M1_beta(fi);
            S_AA_all(:,:,fi) = S_AA_all(:,:,fi) + c * c';
            S_AB_all(:,:,fi) = S_AB_all(:,:,fi) + c * conj(m);
            S_BB_all(:,:,fi) = S_BB_all(:,:,fi) + m * conj(m);
        end

        n_segs    = n_segs + 1;
        start_idx = start_idx + (seg_len - n_overlap);
    end

    % Average over segments and beta frequencies
    S_AA = mean(S_AA_all, 3) / n_segs;
    S_AB = mean(S_AB_all, 3) / n_segs;
    S_BB = mean(S_BB_all, 3) / n_segs;
    S_BB = S_BB(1,1);   % scalar
end

function coh = coh_from_phi(phi, S_AA_reg, S_AA_r, S_AB_r, S_BB)
    S_AB_rot = real(exp(-1i*phi) * S_AB_r);
    w_A      = S_AA_reg \ S_AB_rot;
    num      = abs(w_A' * S_AB_r)^2;
    denom    = real(w_A' * S_AA_r * w_A) * real(S_BB);
    if denom > 0
        coh = num / denom;
    else
        coh = 0;
    end
end

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