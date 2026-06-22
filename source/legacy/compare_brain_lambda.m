%% compare_brain_lambda.m
% Compare brain-EMG DICS coherence with absolute lambda=10 vs lambda='10%'
% for participant 1 (OP00212) only.

clear all; close all; clc;

%% =========================================================================
%  USER CONFIG
%% =========================================================================
fieldtrip_path = 'C:\Users\mspedden\Documents\fieldtrip';
spm_path       = 'C:\Users\mspedden\Documents\spm';
bsc_path       = 'C:\Users\mspedden\Documents\brainspineconnectivity\source';
data_root      = 'C:\spinecoh_data';
save_dir       = 'C:\Users\mspedden\Documents\brainspine_savetest';
geomfile       = 'C:\Leadfields meshes\geometries_cervical_realistic.mat';
T_mat_path     = 'C:\Leadfields meshes\T.mat';

sub            = 'OP00212';
fband          = [10 35];
numpermutation = 500;
fwhm_mm        = 8;
radius_mm      = 3 * (fwhm_mm / 2.355);
rng(1);

lambda_configs = {10, '10%'};
lambda_labels  = {'abs10', 'pct10'};

%% =========================================================================
%  SETUP
%% =========================================================================
addpath(bsc_path); addpath(spm_path);
spm('defaults','EEG');
addpath(fieldtrip_path);
ft_defaults;

fig_dir = fullfile(save_dir, 'figures', 'brain_lambda_compare');
if ~exist(fig_dir,'dir'), mkdir(fig_dir); end

ncol     = 256;
addpath(fullfile(fieldtrip_path,'external','matplotlib'));
cmap_hot = [[0.92 0.92 0.92]; flipud(magma(ncol-1))];

%% =========================================================================
%  LOAD DATA ONCE
%% =========================================================================
fprintf('Loading data for %s...\n', sub);
load(geomfile);   % sources_brain, mesh_brain, mesh_torso etc.

datafile = fullfile(data_root, ['sub-' sub], 'ses-001', 'meg', ...
    'pmergedoe1000mspddfflo45hi45hfcstatic_001_array1.mat');

D       = spm_eeg_load(datafile);
grad_mm = D.sensors('MEG');
ftdat   = spm2fieldtrip(D);

badchans = D.chanlabels(D.badchannels);
cfg = []; cfg.channel = setdiff(ftdat.label, badchans);
ftdat = ft_selectdata(cfg, ftdat);

% Rectify EMG
cfg = []; cfg.rectify = 'yes'; cfg.channel = 'EXG1';
ftdatr = ft_preprocessing(cfg, ftdat);
for k = 1:length(ftdat.trial)
    ftdat.trial{k}(end,:) = ftdatr.trial{k};
end

% Brain sensors only
brainidx  = find(grad_mm.chanpos(:,2) > 200);
braingrad = subset_grad(grad_mm, brainidx);

% Head model and leadfield — computed once, shared across lambda conditions
mesh_brain.unit = 'mm';
cfg = []; cfg.method = 'singleshell';
vol = ft_prepare_headmodel(cfg, mesh_brain);
cfg = []; cfg.sourcemodel = sources_brain; cfg.headmodel = vol;
cfg.grad = braingrad; cfg.reducerank = 'no';
LF = ft_prepare_leadfield(cfg);

% Frequency data — computed once
[statdat, restdat, freqdat] = make_freq_data(ftdat, ftdat.trialinfo, fband);

invpthr = -log10(0.05);
results = struct();

%% =========================================================================
%  LOOP OVER LAMBDA CONDITIONS
%% =========================================================================
for li = 1:numel(lambda_configs)
    lambda = lambda_configs{li};
    label  = lambda_labels{li};
    fprintf('\n--- Lambda: %s ---\n', num2str(lambda));
    t0 = tic;

    %% Common spatial filter — exactly as in run_brain_emg_dics
    cfg_dics = [];
    cfg_dics.grid                  = sources_brain;
    cfg_dics.headmodel             = vol;
    cfg_dics.sourcemodel.leadfield = LF.leadfield;
    cfg_dics.dics.keepfilter       = 'yes';
    cfg_dics.dics.lambda           = lambda;
    cfg_dics.method                = 'dics';
    cfg_dics.refchan               = 'EXG1';
    coh_source = ft_sourceanalysis(cfg_dics, freqdat);


    %% Permutation test — exactly as in run_brain_emg_dics
    cfg_perm = [];
    cfg_perm.grid                  = sources_brain;
    cfg_perm.headmodel             = vol;
    cfg_perm.sourcemodel.leadfield = LF.leadfield;
    cfg_perm.dics.filter           = coh_source.avg.filter;
    cfg_perm.dics.lambda           = lambda;
    cfg_perm.method                = 'dics';
    cfg_perm.refchan               = 'EXG1';
    cfg_perm.permutation           = 'yes';
    cfg_perm.numpermutation        = numpermutation;
source_perm = ft_sourceanalysis(cfg_perm, statdat, restdat);

    %% Build smoother on first iteration
    if li == 1
        fprintf('Building brain smoother (FWHM=%d mm)...\n', fwhm_mm);
        inside_idx    = (1:size(source_perm.pos,1))';
        inside_pos    = source_perm.pos;
        Wsm           = make_gaussian_smoother(inside_pos, fwhm_mm, radius_mm);
        nsourcepoints = size(source_perm.pos,1);
        fprintf('  Inside sources: %d\n', nsourcepoints);
    end

    nPerm = numel(source_perm.trialA);
    [coh_diff, cohDiff_perm] = extract_coh_diff_brain(source_perm, nsourcepoints, nPerm);

    %% Build smoother on first iteration (needs source_perm)
    if li == 1
        fprintf('Building brain smoother (FWHM=%d mm)...\n', fwhm_mm);
        inside_idx    = find(source_perm.inside);
        inside_pos    = source_perm.pos(inside_idx,:);
        Wsm           = make_gaussian_smoother(inside_pos, fwhm_mm, radius_mm);
        nsourcepoints = numel(source_perm.inside);
        fprintf('  Inside sources: %d\n', numel(inside_idx));
    end

    %% Smooth — inside sources only
    coh_diff_inside            = coh_diff(inside_idx);
    cohDiff_perm_inside        = cohDiff_perm(inside_idx,:);
    coh_diff_inside            = Wsm * coh_diff_inside;
    cohDiff_perm_inside        = Wsm * cohDiff_perm_inside;
    coh_diff(inside_idx)       = coh_diff_inside;
    cohDiff_perm(inside_idx,:) = cohDiff_perm_inside;

    %% Threshold on inside sources only
    maxPerm = max(cohDiff_perm_inside, [], 1);
    thr95   = prctile(maxPerm, 95);
    mask    = coh_diff > thr95;

    pvals       = (sum(cohDiff_perm >= coh_diff, 2) + 1) / (nPerm + 1);
    invp        = -log10(pvals);
    invp_masked = invp; invp_masked(~mask) = 0;

    fprintf('  Threshold: %.6f  |  Sig sources: %d / %d\n', ...
        thr95, sum(mask), nsourcepoints);
    fprintf('  Time: %.1f min\n', toc(t0)/60);

    %% MNI peak
    maxval   = max(coh_diff);
    maxidx   = find(coh_diff == maxval);
    maxpos   = source_perm.pos(maxidx,:);
    load(T_mat_path); T_inv = inv(T);
    maxpos_h = [maxpos, ones(size(maxpos,1),1)]';
    x_mni    = (T_inv * maxpos_h)'; x_mni = x_mni(:,1:3);
    fprintf('  Peak MNI: [%.1f %.1f %.1f]\n', x_mni(1,:));

    %% Interpolate onto brain mesh
    source_p     = coh_source;
    source_p.avg.coh = invp_masked;
    cfg_interp = []; cfg_interp.parameter = 'coh';
    cfg_interp.interpmethod = 'sphere_avg';
    cfg_interp.sphereradius = 10;
    brain_int = ft_sourceinterpolate(cfg_interp, source_p, mesh_brain);

    %% Store
    results(li).label     = label;
    results(li).lambda    = lambda;
    results(li).coh_diff  = coh_diff;
    results(li).thr95     = thr95;
    results(li).mask      = mask;
    results(li).invp      = invp;
    results(li).brain_int = brain_int;
    results(li).x_mni     = x_mni;

    save(fullfile(save_dir, sprintf('brainResult_sub%s_%s.mat', sub, label)), ...
        'coh_diff','cohDiff_perm','thr95','mask','invp');

    %% Figure
    hfig = figure;
    cfg_plot = []; cfg_plot.figure = 'gcf'; cfg_plot.method = 'surface';
    cfg_plot.funparameter = 'coh'; cfg_plot.funcolormap = cmap_hot;
    if any(mask)
        cfg_plot.funcolorlim = [invpthr max(invp(mask))];
    else
        cfg_plot.funcolorlim = [invpthr invpthr+1];
    end
    cfg_plot.projmethod = 'nearest'; cfg_plot.surffile = mesh_brain;
    ft_sourceplot(cfg_plot, brain_int);
    view(176,-10); camlight; ax = gca; ax.FontSize = 14;
    hpatch = findobj(gcf,'Type','patch'); set(hpatch,'FaceAlpha',0.9);
    title(sprintf('%s — brain-EMG coherence  lambda=%s', sub, num2str(lambda)), ...
        'Interpreter','none');
    savefig(hfig, fullfile(fig_dir, sprintf('brain_coh_%s_lambda%s.fig', sub, label)));
    saveas(hfig,  fullfile(fig_dir, sprintf('brain_coh_%s_lambda%s.png', sub, label)));
    close(hfig);
end

%% =========================================================================
%  COMPARISON SUMMARY
%% =========================================================================
fprintf('\n=== SUMMARY ===\n');
for li = 1:numel(lambda_configs)
    r = results(li);
    fprintf('  lambda=%-6s  thr=%.4e  sig=%d  peak MNI=[%.1f %.1f %.1f]\n', ...
        num2str(r.lambda), r.thr95, sum(r.mask), r.x_mni(1,:));
end
fprintf('\nDone.\n');

%% =========================================================================
%  LOCAL FUNCTIONS
%% =========================================================================
function [coh_diff, cohDiff_perm] = extract_coh_diff_brain(source_perm, nsourcepoints, nPerm)
    cohDiff_perm = zeros(nsourcepoints, nPerm);
    for i = 1:nPerm
        cohDiff_perm(:,i) = source_perm.trialA(i).coh - source_perm.trialB(i).coh;
    end
    coh_diff = source_perm.avgA.coh - source_perm.avgB.coh;
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

function [statdat, restdat, freqdat] = make_freq_data(dat, trialinfo, fband)
    cfg = []; cfg.output = 'powandcsd'; cfg.method = 'mtmfft';
    cfg.foilim = fband; cfg.tapsmofrq = 1; cfg.keeptrials = 'yes';
    freqdat_tr = ft_freqanalysis(cfg, dat);
    cfg = []; cfg.avgoverfreq = 'yes';
    freqdat_tr = ft_selectdata(cfg, freqdat_tr);

    statidx = find(trialinfo==1);
    restidx = find(trialinfo==2);
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