%% DICS brain-EMG coherence — BEM v2
% Computes coherence between brain sources and EMG reference channel.
% Permutation test (static vs rest) with family-wise multiple comparison
% correction. Optional spatial smoothing of source maps before inference.

clear all
close all
clc

addpath('C:\Users\mspedden\Documents\brainspineconnectivity\source')
addpath('C:\Users\mspedden\Documents\spm')
spm('defaults','EEG')

addpath('C:\Users\mspedden\Documents\fieldtrip')
ft_defaults

save_dir = 'C:\Users\mspedden\Documents\brainspine_save_bemv2';
if ~exist(save_dir,'dir')
    mkdir(save_dir)
end

% OP00220 excluded: very small head (sensors too far)
% OP00226 excluded: could not close headcast
subs = {'OP00212','OP00213','OP00215','OP00219', ...
        'OP00225','OP00221','OP00224'};

generic_dir = 'C:\Users\mspedden\Documents\new_leadfields_and_geom';
geomfile    = fullfile(generic_dir, 'geometries_cervical_realistic.mat');
rng(1)

rectify        = 1;
mult_comp_corr = 1;
fband          = [10 35];

%% ---------------- USER OPTION: spatial smoothing ----------------
doSmooth = 1;     % 0 = no smoothing, 1 = smooth source maps before inference

% brain smoothing applied to source maps
% radius is 3*sigma where sigma = FWHM/2.355
brain_smooth_fwhm_mm   = 8;
brain_smooth_radius_mm = 3 * (brain_smooth_fwhm_mm / 2.355);  % ~10.2 mm

if doSmooth
    result_suffix = sprintf('brainEMG_brainSmooth_%dmm', brain_smooth_fwhm_mm);
else
    result_suffix = 'brainEMG';
end
%% ---------------------------------------------------------------

subjResults = struct();

for ss = 1:length(subs)

    sub = subs{ss};
    fprintf('Processing subject %s (%d/%d)\n', sub, ss, length(subs));

    if strcmp(sub, 'OP00224')
        datwithEMGmerged = fullfile('C:\Users\mspedden\Documents', ...
            ['sub-' sub], 'ses-001', 'meg', ...
            'pmergedoe1000mspddfflo45hi45hfcstatic_002_array1.mat');
    else
        datwithEMGmerged = fullfile('C:\Users\mspedden\Documents', ...
            ['sub-' sub], 'ses-001', 'meg', ...
            'pmergedoe1000mspddfflo45hi45hfcstatic_001_array1.mat');
    end

    load(geomfile)
    D       = spm_eeg_load(datwithEMGmerged);
    grad_mm = D.sensors('MEG');
    ftdat   = spm2fieldtrip(D);

    badchans = D.chanlabels(D.badchannels);
    cfg = []; cfg.channel = setdiff(ftdat.label, badchans);
    ftdat = ft_selectdata(cfg, ftdat);

    if rectify
        cfg = []; cfg.rectify = 'yes'; cfg.channel = 'EXG1';
        ftdatr = ft_preprocessing(cfg, ftdat);
        for k = 1:length(ftdat.trial)
            ftdat.trial{k}(end,:) = ftdatr.trial{k};
        end
    end

    %% brain channel subset
    brainidx           = find(grad_mm.chanpos(:,2) > 200);
    braingrad          = grad_mm;
    braingrad.chanori  = grad_mm.chanori(brainidx, :);
    braingrad.chanpos  = grad_mm.chanpos(brainidx, :);
    braingrad.chantype = grad_mm.chantype(1:length(brainidx));
    braingrad.chanunit = grad_mm.chanunit(1:length(brainidx));
    braingrad.coilori  = grad_mm.coilori(brainidx, :);
    braingrad.coilpos  = grad_mm.coilpos(brainidx, :);
    braingrad.label    = grad_mm.label(brainidx);
    braingrad.tra      = grad_mm.tra(brainidx, brainidx);

    % select brain channels from ftdat and append EMG back for coherence ref
    cfg = []; cfg.channel = braingrad.label;
    braindat = ft_selectdata(cfg, ftdat);
    cfg = []; cfg.channel = {'EXG1'};
    emgdat = ft_selectdata(cfg, ftdat);
    braindat = ft_appenddata([], braindat, emgdat);
    braindat.trialinfo = ftdat.trialinfo;

    mesh_brain.unit = 'mm';

    %% volume conductor + leadfield
    cfg = []; cfg.method = 'singleshell';
    vol = ft_prepare_headmodel(cfg, mesh_brain);

    cfg = []; cfg.sourcemodel = sources_brain; cfg.headmodel = vol;
    cfg.grad = braingrad; cfg.reducerank = 'no';
    LF = ft_prepare_leadfield(cfg);

    %% trial-wise frequency data
    cfg = []; cfg.output = 'powandcsd'; cfg.method = 'mtmfft';
    cfg.foilim = fband; cfg.tapsmofrq = 1; cfg.keeptrials = 'yes';
    freqdat_tr = ft_freqanalysis(cfg, braindat);
    cfg = []; cfg.avgoverfreq = 'yes';
    freqdat_tr = ft_selectdata(cfg, freqdat_tr);

    %% separate conditions
    statidx = find(ftdat.trialinfo==1);
    restidx = find(ftdat.trialinfo==2);
    [nTrials,~] = min([length(statidx) length(restidx)]);

    cfg = []; cfg.trials = statidx(1:nTrials);
    statdat = ft_selectdata(cfg, freqdat_tr);
    cfg = []; cfg.trials = restidx(1:nTrials);
    restdat = ft_selectdata(cfg, freqdat_tr);

    %% mean freq data for filter
    cfg = []; cfg.output = 'powandcsd'; cfg.method = 'mtmfft';
    cfg.foilim = fband; cfg.tapsmofrq = 1; cfg.keeptrials = 'no';
    freqdat = ft_freqanalysis(cfg, braindat);
    cfg = []; cfg.avgoverfreq = 'yes';
    freqdat = ft_selectdata(cfg, freqdat);

    %% DICS filter
    cfg = []; cfg.grid = sources_brain; cfg.headmodel = vol;
    cfg.sourcemodel.leadfield = LF; cfg.dics.keepfilter = 'yes';
    cfg.dics.lambda = 10; cfg.method = 'dics'; cfg.refchan = 'EXG1';
    coh_source = ft_sourceanalysis(cfg, freqdat);

    %% permutation test
    cfg = []; cfg.grid = sources_brain; cfg.headmodel = vol;
    cfg.sourcemodel.leadfield = LF;
    cfg.dics.filter = coh_source.avg.filter;
    cfg.dics.lambda = 10; cfg.method = 'dics'; cfg.refchan = 'EXG1';
    cfg.permutation = 'yes'; cfg.numpermutation = 500;
    source_perm = ft_sourceanalysis(cfg, statdat, restdat);

    nsourcepoints = length(source_perm.inside);
    nPerm         = numel(source_perm.trialA);
    cohDiff_perm  = zeros(nsourcepoints, nPerm);
    for i = 1:nPerm
        cohDiff_perm(:,i) = source_perm.trialA(i).coh - source_perm.trialB(i).coh;
    end

    coh_diff = source_perm.avgA.coh - source_perm.avgB.coh;

    %% optional spatial smoothing of observed + permutation maps
    if doSmooth
        inside_pos = sources_brain(source_perm.inside, :);  % [nInside x 3]
        Wsm = make_gaussian_smoother(inside_pos, brain_smooth_fwhm_mm, brain_smooth_radius_mm);

        nnz_per_row = full(sum(Wsm > 0, 2));
        selfw       = full(diag(Wsm));
        fprintf('  Brain smoother neighbours/row: median %.1f (min %d, max %d)\n', ...
            median(nnz_per_row), min(nnz_per_row), max(nnz_per_row));
        fprintf('  Brain smoother self-weight:    median %.3f (min %.3f, max %.3f)\n', ...
            median(selfw), min(selfw), max(selfw));

        cohDiff_perm_inf = Wsm * cohDiff_perm;
        coh_diff_inf     = Wsm * coh_diff;
    else
        cohDiff_perm_inf = cohDiff_perm;
        coh_diff_inf     = coh_diff;
    end

    maxPerm = max(cohDiff_perm_inf, [], 1);

    if mult_comp_corr
        thr95 = prctile(maxPerm, 95);
    else
        thr95 = prctile(cohDiff_perm_inf, 95, 2);
        warning('Uncorrected per-source threshold used.')
    end

    %% p-values
    pvals = zeros(nsourcepoints,1);
    for s = 1:nsourcepoints
        permDist = cohDiff_perm_inf(s,:);
        obsVal   = coh_diff_inf(s);
        pvals(s) = (sum(permDist >= obsVal) + 1) / (nPerm + 1);
    end

    invp               = -log10(pvals);
    mask               = coh_diff_inf > thr95;
    invp_masked        = invp;
    invp_masked(~mask) = 0;
    invpthr            = -log10(0.05);

    %% MNI coordinates of max
    source_diff         = coh_source;
    source_diff.avg.coh = coh_diff_inf;
    maxval     = max(source_diff.avg.coh);
    maxidx     = find(source_diff.avg.coh == maxval);
    maxpos     = source_diff.pos(maxidx, :);
    load(fullfile(save_dir, 'T.mat'));
    T_inv    = inv(T);
    maxpos_h = [maxpos, ones(size(maxpos,1),1)]';
    x_mni    = (T_inv * maxpos_h)'; x_mni = x_mni(:,1:3);
    if length(maxidx) > 1, disp('multiple max locs'); end

    figure; ft_plot_mesh(mesh_brain); hold on
    plot3(maxpos(:,1), maxpos(:,2), maxpos(:,3), 'r*')
    if doSmooth
        title(sprintf('%s max location (brain smoothed %d mm)', sub, brain_smooth_fwhm_mm), ...
            'Interpreter','none')
    else
        title(sprintf('%s max location', sub), 'Interpreter','none')
    end

    %% interpolate onto brain mesh — sphere_avg for smoother surface map
    source_p         = coh_source;
    source_p.avg.coh = invp_masked;
    cfg              = [];
    cfg.parameter    = 'coh';
    cfg.interpmethod = 'sphere_avg';
    cfg.sphereradius = 10;
    brain_int        = ft_sourceinterpolate(cfg, source_p, mesh_brain);

    source_mask         = coh_source;
    source_mask.avg.coh = double(mask);
    cfg                 = [];
    cfg.parameter       = 'coh';
    cfg.interpmethod    = 'sphere_avg';
    cfg.sphereradius    = 10;
    source_mask_int     = ft_sourceinterpolate(cfg, source_mask, mesh_brain);
    brain_int.mask      = source_mask_int.coh;

    %% colormap
    ncol        = 256;
    addpath('C:\Users\mspedden\Documents\fieldtrip\external\matplotlib\')
    brain_color = [0.92 0.92 0.92];
    hotmap      = flipud(magma(ncol-1));
    cmap        = [brain_color; hotmap];

    %% plot individual
    figure;
    cfg              = [];
    cfg.figure       = 'gcf';
    cfg.method       = 'surface';
    cfg.funparameter = 'coh';
    cfg.funcolormap  = cmap;
    if any(mask)
        cfg.funcolorlim = [invpthr max(invp(mask))];
    else
        cfg.funcolorlim = [invpthr invpthr+1];
    end
    cfg.projmethod = 'nearest';
    cfg.surffile   = mesh_brain;
    ft_sourceplot(cfg, brain_int);
    view(176,-10); camlight;
    ax = gca; ax.FontSize = 14;
    hpatch = findobj(gcf,'Type','patch');
    set(hpatch,'FaceAlpha',0.9)
    if doSmooth
        title(sprintf('%s: brain-EMG coherence (smoothed %d mm)', sub, brain_smooth_fwhm_mm), ...
            'Interpreter','none')
    else
        title(sprintf('%s: brain-EMG coherence', sub), 'Interpreter','none')
    end

    %% store results
    subjResults(ss).coh_diff   = coh_diff_inf;
    subjResults(ss).sig_mask   = source_mask_int.coh;
    subjResults(ss).pos        = source_mask_int.pos;
    subjResults(ss).inside     = source_mask_int.inside;
    subjResults(ss).thr95      = thr95;
    subjResults(ss).maxpos     = maxpos;
    subjResults(ss).maxpos_mni = x_mni;
    subjResults(ss).brain_int  = brain_int;  % save for quick replotting

end

%% save
groupfile = sprintf('groupRes_brain_DICS_bemv2_%s.mat', result_suffix);
save(fullfile(save_dir, groupfile), 'subjResults')

%%------------- GROUP ANALYSIS -------------------

nSubs     = length(subjResults);
all_masks = cat(2, subjResults(:).sig_mask);
sig_pos   = false(nSubs,1);

for ss = 1:nSubs
    if any(subjResults(ss).coh_diff > subjResults(ss).thr95)
        sig_pos(ss) = true;
    end
end

if doSmooth
    fprintf('%g/%g subjects show sig brain-EMG coherence (smoothed %d mm)\n', ...
        sum(sig_pos), nSubs, brain_smooth_fwhm_mm)
else
    fprintf('%g/%g subjects show sig brain-EMG coherence\n', sum(sig_pos), nSubs)
end

out_brain = plot_bayesprev_hpdi_only(sig_pos, 0.05);

x = sum(sig_pos); n = nSubs; alpha = 0.05;
lower = betainv(alpha/2,   x,   n-x+1);
upper = betainv(1-alpha/2, x+1, n-x);
ci = [lower upper];
p  = 1 - binocdf(x-1, n, alpha);
fprintf('Binomial CI: [%.3f %.3f], p = %.4f\n', ci(1), ci(2), p);

%% group prevalence
group_prevalence = mean(all_masks, 2);
threshold        = 0.3;
group_prevalence_masked = group_prevalence;
group_prevalence_masked(group_prevalence < threshold) = 0;

group_ft        = [];
group_ft.pos    = subjResults(1).pos;
group_ft.inside = subjResults(1).inside;
group_ft.pow    = group_prevalence_masked;

%% MNI max for group
maxval   = max(group_ft.pow);
maxidx   = find(group_ft.pow == maxval);
maxpos   = group_ft.pos(maxidx, :);
maxpos_h = [maxpos, ones(size(maxpos,1),1)]';
x_mni    = (T_inv * maxpos_h)'; x_mni = x_mni(:,1:3);
if length(maxidx) > 1, disp('multiple max locs'); end

figure; ft_plot_mesh(mesh_brain); hold on
plot3(maxpos(:,1), maxpos(:,2), maxpos(:,3), 'r*')
title('Group max location','Interpreter','none')

%% interpolate group prevalence
cfg              = [];
cfg.parameter    = 'pow';
cfg.interpmethod = 'sphere_avg';
cfg.sphereradius = 10;
group_int = ft_sourceinterpolate(cfg, group_ft, mesh_brain);

%% plot group prevalence
figure;
cfg              = [];
cfg.method       = 'surface';
cfg.funparameter = 'pow';
cfg.funcolorlim  = [threshold max(group_int.pow)];
cfg.funcolormap  = cmap;
cfg.projmethod   = 'nearest';
cfg.surffile     = mesh_brain;
ft_sourceplot(cfg, group_int);
view(176,-10); camlight;
ax = gca; ax.FontSize = 14;
hpatch = findobj(gcf,'Type','patch');
set(hpatch,'FaceAlpha',0.9)
if doSmooth
    title(sprintf('Group prevalence — brain-EMG coherence (smoothed %d mm)', ...
        brain_smooth_fwhm_mm), 'Interpreter','none')
else
    title('Group prevalence — brain-EMG coherence', 'Interpreter','none')
end


%% =========================================================
%  HELPER FUNCTION
%  =========================================================

function W = make_gaussian_smoother(pos_mm, fwhm_mm, radius_mm)
% Returns an [N x N] row-normalised sparse smoothing matrix W such that:
%   x_smoothed = W * x
%
% INPUTS:
%   pos_mm    - [N x 3] source positions in mm (inside sources only)
%   fwhm_mm   - smoothing kernel FWHM in mm
%   radius_mm - neighbourhood search radius; should be 3*(fwhm_mm/2.355)

sigma = fwhm_mm / 2.355;

if nargin < 3 || isempty(radius_mm)
    radius_mm = 3 * sigma;
end

N   = size(pos_mm, 1);
Mdl = KDTreeSearcher(pos_mm);
[idx, dist] = rangesearch(Mdl, pos_mm, radius_mm);

ii = [];
jj = [];
vv = [];

for i = 1:N
    j = idx{i};
    d = dist{i};
    w = exp(-0.5 * (d ./ sigma).^2);
    ii = [ii; repmat(i, numel(j), 1)];
    jj = [jj; j(:)];
    vv = [vv; w(:)];
end

W  = sparse(ii, jj, vv, N, N);
rs = full(sum(W, 2));
rs(rs == 0) = 1;
W  = spdiags(1./rs, 0, N, N) * W;

end
