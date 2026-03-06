%% DICS spinal VE to brain — BEM v2
% Functional spinal VE only
%
% Optional spatial smoothing:
% If doSmooth = 1, then:
%   1) load the spinal VE derived from the smoothed spinal analysis
%   2) spatially smooth the brain coherence difference and permutation maps
%   3) save outputs with filenames that encode both:
%        - spinal smoothing used upstream for VE definition
%        - brain smoothing used here for brain inference
%
% Current convention:
%   spinal smoothing (upstream VE source) = 40 mm
%   brain smoothing (this script)         = 20 mm

clear all
close all
clc

addpath('C:\Users\mspedden\Documents\brainspineconnectivity\source')
addpath('C:\Users\mspedden\Documents\spm')
spm('defaults','EEG')

addpath('C:\Users\mspedden\Documents\fieldtrip')
ft_defaults

save_dir='C:\Users\mspedden\Documents\brainspine_save_bemv2';
if ~exist(save_dir,'dir')
    mkdir(save_dir)
end

subs = {'OP00212','OP00213','OP00215', 'OP00219', ...
    'OP00225', 'OP00221', 'OP00224'};

generic_dir = 'C:\Users\mspedden\Documents\new_leadfields_and_geom';
geomfile = fullfile(generic_dir, 'geometries_cervical_realistic.mat');
rng(1)

HFC            = 1;
rectify        = 1;
mult_comp_corr = 1;
fband          = [10 35];

%% ---------------- USER OPTION: spatial smoothing ----------------
doSmooth = 1;     % 0 = no smoothing, 1 = load smoothed spinal VE + smooth brain maps

% spinal smoothing used earlier in pipeline (used only for selecting VE file)
spine_smooth_fwhm_mm = 40;

% brain smoothing applied here to brain source maps
brain_smooth_fwhm_mm   = 20;
brain_smooth_radius_mm = 3*brain_smooth_fwhm_mm;

if doSmooth
    result_suffix = sprintf('functionalVE_spineSmooth_%dmm_brainSmooth_%dmm', ...
        spine_smooth_fwhm_mm, brain_smooth_fwhm_mm);
    loadfile_suffix = sprintf('_permSmooth_%dmm', spine_smooth_fwhm_mm);
else
    result_suffix = 'functionalVE';
    loadfile_suffix = '';
end
%% ---------------------------------------------------------------

subjResults = struct();

for ss = 1:length(subs)

    sub = subs{ss};

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
    cfg = [];
    cfg.channel = setdiff(ftdat.label, badchans);
    ftdat = ft_selectdata(cfg, ftdat);

    if rectify
        cfg = [];
        cfg.rectify = 'yes';
        cfg.channel = 'EXG1';
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

    %% load spinal VE
    loadfile = sprintf('VE_spine_sub%s_forspectra_bemv2%s.mat', sub, loadfile_suffix);
    load(fullfile(save_dir, loadfile), 'VE')

    trialinfo = ftdat.trialinfo;
    braindat  = ft_appenddata([], ftdat, VE);
    braindat.trialinfo = trialinfo;

    mesh_brain.unit = 'mm';
    braindat.grad   = braingrad;

    %% single shell volume conductor
    cfg        = [];
    cfg.method = 'singleshell';
    vol        = ft_prepare_headmodel(cfg, mesh_brain);

    %% brain leadfield
    cfg             = [];
    cfg.sourcemodel = sources_brain;
    cfg.headmodel   = vol;
    cfg.grad        = braingrad;
    cfg.reducerank  = 'no';
    cfg.normalize   = 'no';
    LF = ft_prepare_leadfield(cfg);

    %% trial-wise frequency data
    cfg            = [];
    cfg.output     = 'powandcsd';
    cfg.method     = 'mtmfft';
    cfg.foilim     = fband;
    cfg.tapsmofrq  = 1;
    cfg.keeptrials = 'yes';
    freqdat_tr = ft_freqanalysis(cfg, braindat);

    cfg = [];
    cfg.avgoverfreq = 'yes';
    freqdat_tr = ft_selectdata(cfg, freqdat_tr);

    %% separate conditions
    statidx = find(ftdat.trialinfo==1);
    restidx = find(ftdat.trialinfo==2);
    [nTrials,~] = min([length(statidx) length(restidx)]);

    cfg = [];
    cfg.trials = statidx(1:nTrials);
    statdat = ft_selectdata(cfg, freqdat_tr);

    cfg = [];
    cfg.trials = restidx(1:nTrials);
    restdat = ft_selectdata(cfg, freqdat_tr);

    %% mean freq data for filter
    cfg            = [];
    cfg.output     = 'powandcsd';
    cfg.method     = 'mtmfft';
    cfg.foilim     = fband;
    cfg.tapsmofrq  = 1;
    cfg.keeptrials = 'no';
    freqdat = ft_freqanalysis(cfg, braindat);

    cfg = [];
    cfg.avgoverfreq = 'yes';
    freqdat = ft_selectdata(cfg, freqdat);

    %% DICS filter
    cfg                       = [];
    cfg.grid                  = sources_brain;
    cfg.headmodel             = vol;
    cfg.sourcemodel.leadfield = LF;
    cfg.dics.keepfilter       = 'yes';
    cfg.dics.lambda           = 10;
    cfg.method                = 'dics';
    cfg.refchan               = 'virtualchannel001';
    coh_source = ft_sourceanalysis(cfg, freqdat);

    %% permutation test
    cfg                       = [];
    cfg.grid                  = sources_brain;
    cfg.headmodel             = vol;
    cfg.sourcemodel.leadfield = LF;
    cfg.dics.filter           = coh_source.avg.filter;
    cfg.dics.lambda           = 10;
    cfg.method                = 'dics';
    cfg.refchan               = 'virtualchannel001';
    cfg.permutation           = 'yes';
    cfg.numpermutation        = 500;
    source_perm = ft_sourceanalysis(cfg, statdat, restdat);

    nsourcepoints = length(source_perm.inside);

    nPerm        = numel(source_perm.trialA);
    cohDiff_perm = zeros(nsourcepoints, nPerm);
    for i = 1:nPerm
        cohDiff_perm(:,i) = source_perm.trialA(i).coh - source_perm.trialB(i).coh;
    end

    coh_diff = source_perm.avgA.coh - source_perm.avgB.coh;

    %% optional spatial smoothing of observed + permutation maps in brain source space
    if doSmooth
        Wsm = make_gaussian_smoother(sources_brain.pos, brain_smooth_fwhm_mm, brain_smooth_radius_mm);

        nnz_per_row = full(sum(Wsm>0,2));
        selfw = full(diag(Wsm));
        fprintf('Brain smoother neighbours/row: median %.1f (min %d, max %d)\n', ...
            median(nnz_per_row), min(nnz_per_row), max(nnz_per_row));
        fprintf('Brain smoother self-weight: median %.3f (min %.3f, max %.3f)\n', ...
            median(selfw), min(selfw), max(selfw));

        cohDiff_perm_inf = Wsm * cohDiff_perm;
        coh_diff_inf     = Wsm * coh_diff;
    else
        cohDiff_perm_inf = cohDiff_perm;
        coh_diff_inf     = coh_diff;
    end

    maxPerm = max(cohDiff_perm_inf, [], 1);

    if mult_comp_corr
        thr95 = prctile(maxPerm, 95, 2);
    else
        thr95 = prctile(cohDiff_perm_inf, 95, 2);
    end

    %% p-values
    pvals = zeros(nsourcepoints,1);
    for s = 1:nsourcepoints
        permDist = cohDiff_perm_inf(s,:);
        obsVal   = coh_diff_inf(s);
        pvals(s) = (sum(permDist >= obsVal) + 1) / (nPerm + 1);
    end

    invp    = -log10(pvals);
    mask    = coh_diff_inf > thr95;
    invp_masked        = invp;
    invp_masked(~mask) = 0;
    pthr    = 0.05;
    invpthr = -log10(pthr);

    %% interpolate onto brain mesh
    source_p         = coh_source;
    source_p.avg.coh = invp_masked;

    cfg = [];
    cfg.parameter = 'coh';
    brain_int = ft_sourceinterpolate(cfg, source_p, mesh_brain);

    source_mask         = coh_source;
    source_mask.avg.coh = double(mask);
    cfg                 = [];
    cfg.parameter       = 'coh';
    cfg.interpmethod    = 'nearest';
    source_mask_int     = ft_sourceinterpolate(cfg, source_mask, mesh_brain);
    brain_int.mask      = source_mask_int.coh;

    %% colormap
    ncol = 256;
    addpath('C:\Users\mspedden\Documents\fieldtrip\external\matplotlib\')
    brain_color = [0.92 0.92 0.92];
    hotmap      = flipud(magma(ncol-1));
    cmap        = [brain_color; hotmap];

    %% plot individual
    figure;
    cfg = [];
    cfg.figure       = 'gcf';
    cfg.method       = 'surface';
    cfg.funparameter = 'coh';
    cfg.funcolormap  = cmap;

    if any(mask)
        cfg.funcolorlim = [invpthr max(invp(mask))];
    else
        cfg.funcolorlim = [invpthr invpthr+1];
    end

    cfg.projmethod   = 'nearest';
    cfg.surffile     = mesh_brain;
    ft_sourceplot(cfg, brain_int);
    view(176,-10);
    camlight;
    ax = gca;
    ax.FontSize = 14;
    hpatch = findobj(gcf, 'Type', 'patch');
    set(hpatch, 'FaceAlpha', 0.9)

    if doSmooth
        title(sprintf('%s: brain coherence with functional spinal VE (brain smoothed %d mm)', ...
            sub, brain_smooth_fwhm_mm), 'Interpreter', 'none')
    else
        title(sprintf('%s: brain coherence with functional spinal VE', ...
            sub), 'Interpreter', 'none')
    end

    %% max location in native + MNI
    source_diff         = coh_source;
    source_diff.avg.coh = coh_diff_inf;

    maxval     = max(source_diff.avg.coh);
    maxcohindx = find(source_diff.avg.coh == maxval);
    maxpos     = source_diff.pos(maxcohindx, :);

    load('C:\Users\mspedden\Documents\brainspine_save\T.mat');
    T_inv    = inv(T);
    maxpos_h = [maxpos, ones(size(maxpos,1),1)]';
    x_mni_h  = T_inv * maxpos_h;
    x_mni    = x_mni_h(1:3,:)';

    if length(maxcohindx) > 1
        disp('multiple max locs')
    end

    figure;
    ft_plot_mesh(mesh_brain);
    hold on
    plot3(maxpos(:,1), maxpos(:,2), maxpos(:,3), 'r*')

    if doSmooth
        title(sprintf('%s max native-space location (brain smoothed %d mm)', ...
            sub, brain_smooth_fwhm_mm), 'Interpreter', 'none')
    else
        title(sprintf('%s max native-space location', sub), 'Interpreter', 'none')
    end

    %% save subject results
    subjResults(ss).subjID      = sub;
    subjResults(ss).coh_diff    = coh_diff_inf;
    subjResults(ss).sig_mask    = source_mask_int.coh;
    subjResults(ss).pos         = brain_int.pos;
    subjResults(ss).inside      = brain_int.inside;
    subjResults(ss).maxpos      = maxpos;
    subjResults(ss).maxpos_mni  = x_mni;
    subjResults(ss).thr95       = thr95;

end

%% save group results
groupfile = sprintf('groupRes_brain_DICS_spineVC_bemv2_%s.mat', result_suffix);
save(fullfile(save_dir, groupfile), 'subjResults')

%%------------- GROUP ANALYSIS -------------------

nSubjects = numel(subjResults);
sig_res   = false(nSubjects,1);

for i = 1:nSubjects
    sig_res(i) = any(subjResults(i).sig_mask(:));
end

if doSmooth
    fprintf('%g/%g subjects show sig coherence anywhere (spine smoothed %d mm, brain smoothed %d mm)\n', ...
        sum(sig_res), nSubjects, spine_smooth_fwhm_mm, brain_smooth_fwhm_mm)
else
    fprintf('%g/%g subjects show sig coherence anywhere\n', sum(sig_res), nSubjects)
end

%% group prevalence
all_masks        = cat(2, subjResults(:).sig_mask);
group_prevalence = mean(all_masks, 2);

group_ft        = [];
group_ft.pos    = subjResults(1).pos;
group_ft.inside = subjResults(1).inside;

threshold = 0.2;
group_prevalence_masked = group_prevalence;
group_prevalence_masked(group_prevalence < threshold) = 0;
group_ft.pow = group_prevalence_masked;

%% mask for volumetric data
mask_vol = false(size(group_ft.pow));
mask_vol(group_ft.inside) = group_ft.pow(group_ft.inside) >= threshold;

%% connected-component clustering
inside_idx    = find(group_ft.inside);
voxel_pos     = group_ft.pos(inside_idx, :);
neighbourdist = 20;

active_idx = find(mask_vol(group_ft.inside));
active_pos = voxel_pos(active_idx, :);
active_N   = numel(active_idx);
labels     = zeros(active_N,1);
current_label = 0;

for i = 1:active_N
    if labels(i) == 0
        current_label = current_label + 1;
        stack     = i;
        labels(i) = current_label;
        while ~isempty(stack)
            v = stack(end);
            stack(end) = [];
            neighbors_v = find(vecnorm(active_pos - active_pos(v,:), 2, 2) <= neighbourdist & ...
                               vecnorm(active_pos - active_pos(v,:), 2, 2) > 0);
            neighbors_v = neighbors_v(labels(neighbors_v) == 0);
            labels(neighbors_v) = current_label;
            stack = [stack; neighbors_v];
        end
    end
end

%% interpolate group prevalence onto brain mesh
cfg = [];
cfg.parameter    = 'pow';
cfg.interpmethod = 'nearest';
group_int = ft_sourceinterpolate(cfg, group_ft, mesh_brain);

%% plot group prevalence
figure;
cfg = [];
cfg.method       = 'surface';
cfg.funparameter = 'pow';
cfg.funcolorlim  = [threshold 0.5];
cfg.funcolormap  = cmap;
cfg.projmethod   = 'nearest';
cfg.surffile     = mesh_brain;
ft_sourceplot(cfg, group_int);
view(176, -10);
camlight('headlight');
hold on;
ax = gca;
ax.FontSize = 14;
hpatch = findobj(gcf, 'Type', 'patch');
set(hpatch, 'FaceAlpha', 0.9)

if doSmooth
    title(sprintf('Group prevalence (functional VE; spine %d mm, brain %d mm)', ...
        spine_smooth_fwhm_mm, brain_smooth_fwhm_mm), 'Interpreter', 'none')
else
    title('Group prevalence (functional VE)', 'Interpreter', 'none')
end

%% max location in MNI
maxval   = max(group_ft.pow);
maxidx   = find(group_ft.pow == maxval);
maxpos   = group_ft.pos(maxidx, :);
maxpos_h = [maxpos, ones(size(maxpos,1),1)]';
x_mni_h  = T_inv * maxpos_h;
x_mni    = x_mni_h(1:3,:)';

figure;
ft_plot_mesh(mesh_brain);
hold on
plot3(maxpos(:,1), maxpos(:,2), maxpos(:,3), 'r*')

if doSmooth
    title(sprintf('Group max location (brain smoothed %d mm)', brain_smooth_fwhm_mm), ...
        'Interpreter', 'none')
else
    title('Group max location', 'Interpreter', 'none')
end

%% ROI overlap with active voxels
load('C:\Users\mspedden\Documents\brainspine_save\roi_native.mat')

dist_tol         = 1e-3;
D                = pdist2(roi_native, active_pos);
overlap_logical  = D <= dist_tol;
roi_overlap      = any(overlap_logical, 2);
intersection_pos = roi_native(roi_overlap, :);

figure;
ft_plot_mesh(mesh_brain, 'facecolor', [0.8 0.8 0.8], 'facealpha', 0.4, 'edgecolor', 'none');
hold on;
plot3(intersection_pos(:,1), intersection_pos(:,2), intersection_pos(:,3), ...
    'go', 'MarkerSize', 6, 'MarkerFaceColor', 'g', 'LineWidth', 2);
axis equal;
view(90,0);
camlight;
lighting gouraud;

if doSmooth
    title(sprintf('ROI overlap (brain smoothed %d mm)', brain_smooth_fwhm_mm), ...
        'Interpreter', 'none')
else
    title('ROI overlap', 'Interpreter', 'none')
end

%% save ROI
roifile = sprintf('M1_ROI_bemv2_%s.mat', result_suffix);
save(fullfile(save_dir, roifile), 'intersection_pos')

%% helper: Gaussian spatial smoothing operator
function W = make_gaussian_smoother(pos_mm, fwhm_mm, radius_mm)
% Returns an [N x N] row-normalised sparse smoothing matrix W such that:
%   x_sm = W * x

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

rs = full(sum(W,2));
rs(rs==0) = 1;
W = spdiags(1./rs, 0, N, N) * W;
end