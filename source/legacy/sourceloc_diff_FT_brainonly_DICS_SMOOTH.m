%% script that runs DICS in FT to image brain coherence with ref chan (EMG)
% + spatial smoothing applied to OBSERVED + PERMUTATIONS (as in spine script)
%
% Key changes:
%  1) New save_dir (won't overwrite)
%  2) Build Gaussian smoother Wsm on INSIDE grid points (sources_brain positions)
%  3) Smooth coh_diff and cohDiff_perm BEFORE maxPerm/thr95, mask, p-values
%  4) Store smoothed coh_diff in subjResults (what inference used)

clear all
close all
clc

addpath('C:\Users\mspedden\Documents\brainspineconnectivity\source')
addpath('C:\Users\mspedden\Documents\spm')
spm('defaults','EEG')

addpath('C:\Users\mspedden\Documents\fieldtrip')
ft_defaults

% -------------------- NEW: versioned output folder -------------------------
base_save_dir = 'C:\Users\mspedden\Documents\brainspine_save';
smooth_fwhm_mm   = 20;                 % brain smoothing
smooth_radius_mm = 3*smooth_fwhm_mm;    % robust support

save_dir = fullfile(base_save_dir, sprintf('brain_EMG_permSmooth_%dmm', smooth_fwhm_mm));
if ~exist(save_dir,'dir'); mkdir(save_dir); end
% --------------------------------------------------------------------------

subs = {'OP00212','OP00213','OP00215', 'OP00219', ...
        'OP00225', 'OP00221', 'OP00224'};

generic_dir = 'C:\Users\mspedden\Documents\new_leadfields_and_geom';
geomfile = fullfile(generic_dir, 'geometries_cervical_realistic.mat');

rng(1)

HFC=1;
rectify=1;
mult_comp_corr=1;
fband=[10 35];

subjResults=struct();

for ss=1:length(subs)

    sub=subs{ss};

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
    D = spm_eeg_load(datwithEMGmerged);
    grad_mm = D.sensors('MEG');
    ftdat = spm2fieldtrip(D);

    badchans = D.chanlabels(D.badchannels);

    % remove bad channels
    cfg=[];
    cfg.channel=setdiff(ftdat.label,badchans);
    ftdat=ft_selectdata(cfg,ftdat);

    % rectify EMG
    if rectify
        cfg=[];
        cfg.rectify='yes';
        cfg.channel='EXG1';
        ftdatr=ft_preprocessing(cfg,ftdat);

        for k=1:length(ftdat.trial)
            ftdat.trial{k}(end,:)=ftdatr.trial{k};
        end
    end

    %% brain-only grad
    brainidx = find(grad_mm.chanpos(:,2) > 200);

    braingrad = grad_mm;
    braingrad.chanori  = grad_mm.chanori(brainidx, :);
    braingrad.chanpos  = grad_mm.chanpos(brainidx, :);
    braingrad.chantype = grad_mm.chantype(brainidx);
    braingrad.chanunit = grad_mm.chanunit(brainidx);
    braingrad.coilori  = grad_mm.coilori(brainidx, :);
    braingrad.coilpos  = grad_mm.coilpos(brainidx, :);
    braingrad.label    = grad_mm.label(brainidx);
    braingrad.tra      = grad_mm.tra(brainidx, brainidx);

    % Use all trials but brain grad
    braindat = ftdat;
    mesh_brain.unit='mm';
    braindat.grad = braingrad;

    %% headmodel (single shell)
    cfg = [];
    cfg.method = 'singleshell';
    vol = ft_prepare_headmodel(cfg, mesh_brain);

    %% leadfield
    cfg                     = [];
    cfg.sourcemodel         = sources_brain;
    cfg.headmodel           = vol;
    cfg.grad                = braingrad;
    cfg.reducerank          = 'no';
    cfg.normalize           = 'no';
    LF = ft_prepare_leadfield(cfg);

    %% trialwise freq
    cfg=[];
    cfg.output     = 'powandcsd';
    cfg.method     = 'mtmfft';
    cfg.foilim     = fband;
    cfg.tapsmofrq  = 1;
    cfg.keeptrials = 'yes';
    freqdat_tr = ft_freqanalysis(cfg, braindat);

    cfg=[];
    cfg.avgoverfreq='yes';
    freqdat_tr = ft_selectdata(cfg, freqdat_tr);

    %% separate conditions
    statidx=find(ftdat.trialinfo==1);
    restidx=find(ftdat.trialinfo==2);

    [nTrials,~] = min([length(statidx) length(restidx)]);

    cfg=[];
    cfg.trials=statidx(1:nTrials);
    statdat=ft_selectdata(cfg,freqdat_tr);

    cfg.trials=restidx(1:nTrials);
    restdat=ft_selectdata(cfg,freqdat_tr);

    %% mean freq (for filter)
    cfg=[];
    cfg.output     = 'powandcsd';
    cfg.method     = 'mtmfft';
    cfg.foilim     = fband;
    cfg.tapsmofrq  = 1;
    cfg.keeptrials = 'no';
    % NOTE: original script used ftdat here; keep braindat for consistency
    freqdat = ft_freqanalysis(cfg, braindat);

    cfg=[];
    cfg.avgoverfreq='yes';
    freqdat=ft_selectdata(cfg,freqdat);

    %% DICS filters on combined data
    cfg=[];
    cfg.grid = sources_brain;
    cfg.headmodel=vol;
    cfg.sourcemodel.leadfield=LF;
    cfg.dics.keepfilter='yes';
    cfg.dics.lambda=10;
    cfg.method = 'dics';
    cfg.refchan='EXG1';

    coh_source = ft_sourceanalysis(cfg,freqdat);

    %% Apply filters + permutation
    cfg=[];
    cfg.grid = sources_brain;
    cfg.headmodel=vol;
    cfg.sourcemodel.leadfield=LF;
    cfg.dics.filter=coh_source.avg.filter;
    cfg.dics.lambda=10;
    cfg.method = 'dics';
    cfg.refchan='EXG1';

    cfg.permutation = 'yes';
    cfg.numpermutation=500;

    rng(1); % reproducible permutations regardless of other RNG usage
    source_perm = ft_sourceanalysis(cfg, statdat, restdat);

    % -------------------- INSIDE indexing -------------------------
    inside = source_perm.inside(:);
    inside_idx = find(inside);
    nInside = numel(inside_idx);

    coh_diff_all = source_perm.avgA.coh - source_perm.avgB.coh;
    coh_diff_in  = coh_diff_all(inside_idx);

    nPerm = numel(source_perm.trialA);
    cohDiff_perm_in = zeros(nInside, nPerm);

    for p = 1:nPerm
        d_all = source_perm.trialA(p).coh - source_perm.trialB(p).coh;
        cohDiff_perm_in(:,p) = d_all(inside_idx);
    end
    % -------------------------------------------------------------

    % -------------------- NEW: spatial smoothing ------------------
    pos_in = source_perm.pos(inside_idx,:); % mm
    Wsm = make_gaussian_smoother(pos_in, smooth_fwhm_mm, smooth_radius_mm);

    cohDiff_perm_sm_in = Wsm * cohDiff_perm_in;
    coh_diff_sm_in     = Wsm * coh_diff_in;
    % -------------------------------------------------------------

    % Max-stat threshold from SMOOTHED permutations (inside only)
    maxPerm = max(cohDiff_perm_sm_in, [], 1);

    if mult_comp_corr
        thr95 = prctile(maxPerm, 95);   % scalar FWE threshold
    else
        thr95 = prctile(cohDiff_perm_sm_in, 95, 2); % per-voxel (uncorrected)
        warning('Uncorrected threshold used (per-voxel).')
    end

    % Mask on SMOOTHED observed
    mask_in = coh_diff_sm_in > thr95;

    % One-sided p-values from SMOOTHED distributions
    pvals_in = zeros(nInside,1);
    for s = 1:nInside
        permDist = cohDiff_perm_sm_in(s,:);
        obsVal   = coh_diff_sm_in(s);
        pvals_in(s) = (sum(permDist >= obsVal) + 1) / (nPerm + 1);
    end

    invp_in = -log10(pvals_in);

    invp_masked_in = invp_in;
    invp_masked_in(~mask_in) = 0;

    pthr = 0.05;
    invpthr = -log10(pthr);

    % Reinsert into ALL-grid vectors for interpolation/plotting
    invp_masked_all = zeros(size(source_perm.pos,1),1);
    invp_masked_all(inside_idx) = invp_masked_in;

    mask_all = zeros(size(source_perm.pos,1),1);
    mask_all(inside_idx) = double(mask_in);

    coh_diff_sm_all = zeros(size(source_perm.pos,1),1);
    coh_diff_sm_all(inside_idx) = coh_diff_sm_in;

    %% Max location (use SMOOTHED coh-diff used for inference)
    if any(mask_in)
        tmp = coh_diff_sm_in; tmp(~mask_in) = -inf;
        [~, imax_in] = max(tmp);
    else
        [~, imax_in] = max(coh_diff_sm_in);
    end
    maxidx_all = inside_idx(imax_in);
    maxpos = source_perm.pos(maxidx_all,:);

    % Report MNI for max location
    load(fullfile(base_save_dir, 'T.mat'));
    T_inv=inv(T);
    maxpos_h = [maxpos, 1]';
    x_mni_h = T_inv * maxpos_h;
    x_mni = x_mni_h(1:3,:)';

    if numel(maxidx_all)>1
        disp('multiple max locs')
    end

    figure; ft_plot_mesh(mesh_brain);
    hold on
    plot3(maxpos(:,1), maxpos(:,2), maxpos(:,3), 'r*')
    title(sprintf('Max (smoothed) native [%0.1f %0.1f %0.1f], MNI [%0.1f %0.1f %0.1f]', ...
        maxpos(1),maxpos(2),maxpos(3), x_mni(1),x_mni(2),x_mni(3)));

    %% Build source struct for p-values
    source_p = coh_source;
    source_p.avg.coh = invp_masked_all;

    cfg = [];
    cfg.parameter = 'coh';
    brain_int = ft_sourceinterpolate(cfg, source_p, mesh_brain);

    %% Interpolated mask
    source_mask = coh_source;
    source_mask.avg.coh = mask_all;

    cfg = [];
    cfg.parameter = 'coh';
    cfg.interpmethod = 'nearest';
    source_mask_int = ft_sourceinterpolate(cfg, source_mask, mesh_brain);

    % Optional: attach mask to brain_int if you want to use cfg.maskparameter
    % brain_int.mask = source_mask_int.coh;

    %% Colormap
    ncol = 256;
    addpath('C:\Users\mspedden\Documents\fieldtrip\external\matplotlib\')
    brain_color = [0.92 0.92 0.92];
    hotmap = flipud(magma(ncol-1));
    cmap = [brain_color; hotmap];

    %% Plot -log10(p) masked
    figure;
    cfg = [];
    cfg.figure='gcf';
    cfg.method = 'surface';
    cfg.funparameter = 'coh';
    cfg.funcolormap = cmap;

    sigvals = invp_in(mask_in & isfinite(invp_in));
    if isempty(sigvals)
        cfg.funcolorlim = [invpthr invpthr+1];
    else
        cfg.funcolorlim = [invpthr max(sigvals)];
    end

    cfg.projmethod = 'nearest';
    cfg.surffile = mesh_brain;   % FIXED: use variable, not string
    ft_sourceplot(cfg, brain_int);
    view(176,-10)
    camlight
    ax = gca;
    ax.FontSize = 14;
    hpatch = findobj(gcf, 'Type', 'patch');
    set(hpatch, 'FaceAlpha',0.9)

    %% Save subject results (SMOOTHED inference)
    subjResults(ss).coh_diff = coh_diff_sm_all;     % stored on full grid
    subjResults(ss).sig_mask = source_mask_int.coh; % interpolated binary (0/1) on mesh
    subjResults(ss).pos      = source_mask_int.pos;
    subjResults(ss).inside   = source_mask_int.inside;
    subjResults(ss).thr95    = thr95;
    subjResults(ss).subjID   = sub;

    subjResults(ss).maxdiff.idx_all = maxidx_all;
    subjResults(ss).maxdiff.pos_native = maxpos;
    subjResults(ss).maxdiff.pos_mni = x_mni;

end

save(fullfile(save_dir,'groupRes_brain_DICS_permSmooth.mat'), 'subjResults')

%% ---------------- GROUP / PREVALENCE ----------------
load(fullfile(save_dir,'groupRes_brain_DICS_permSmooth.mat'), 'subjResults')

nSubs = length(subjResults);

% cat all binary significance masks (already interpolated onto mesh)
all_masks = cat(2, subjResults(:).sig_mask);

sig_pos = false(nSubs,1);
for ss = 1:nSubs
    % significant if any vertex is significant
    sig_pos(ss) = any(subjResults(ss).sig_mask(:) > 0);
end

fprintf('%g out of %g subjects show sig coherence in brain (smoothed inference)\n', sum(sig_pos), nSubs)

% out_brain = plot_bayesprev_hpdi_only(sig_pos, 0.05);

%% prevalence (on mesh interpolation grid)
group_prevalence = mean(all_masks, 2);

group_ft = [];
group_ft.pos    = subjResults(1).pos;
group_ft.inside = subjResults(1).inside;
group_ft.pow    = group_prevalence;

threshold = 0.3;
group_ft.pow(group_ft.pow < threshold) = 0;

cfg = [];
cfg.parameter = 'pow';
cfg.interpmethod = 'nearest';
group_int = ft_sourceinterpolate(cfg, group_ft, mesh_brain);

%% Plot group prevalence map
figure;
cfg = [];
cfg.method       = 'surface';
cfg.funparameter = 'pow';
cfg.funcolorlim  = [threshold max(group_int.pow)];
cfg.funcolormap  = cmap;
cfg.projmethod   = 'nearest';
cfg.surffile     = mesh_brain;
ft_sourceplot(cfg, group_int);
view(176, -10);
camlight;
ax = gca;
ax.FontSize = 14;
hpatch = findobj(gcf, 'Type', 'patch');
set(hpatch, 'FaceAlpha',0.9)

%% ------------------------------------------------------------------------
% Helper: Gaussian spatial smoothing operator (sparse, row-normalised)
function W = make_gaussian_smoother(pos_mm, fwhm_mm, radius_mm)
% Returns an [N x N] row-normalised sparse smoothing matrix W such that:
%   x_sm = W * x
%
% pos_mm: [N x 3] positions in mm
% fwhm_mm: scalar FWHM in mm
% radius_mm: truncate kernel beyond this radius

if nargin < 3 || isempty(radius_mm)
    sigma = fwhm_mm/2.355;
    radius_mm = 3*sigma;
else
    sigma = fwhm_mm/2.355;
end

N = size(pos_mm,1);

Mdl = KDTreeSearcher(pos_mm);
[idx, dist] = rangesearch(Mdl, pos_mm, radius_mm);

ii = [];
jj = [];
vv = [];

for i = 1:N
    j = idx{i};
    d = dist{i};
    w = exp(-0.5*(d./sigma).^2);

    ii = [ii; repmat(i, numel(j), 1)];
    jj = [jj; j(:)];
    vv = [vv; w(:)];
end

W = sparse(ii, jj, vv, N, N);

% Row-normalise
rs = full(sum(W,2));
rs(rs==0) = 1;
W = spdiags(1./rs, 0, N, N) * W;
end