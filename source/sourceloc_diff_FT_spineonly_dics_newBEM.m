% SPINE EMG COH

clear all
close all
clc

addpath('C:\Users\mspedden\Documents\brainspineconnectivity\source')
addpath('C:\Users\mspedden\Documents\spm')
spm('defaults','EEG')

addpath('C:\Users\mspedden\Documents\fieldtrip')
ft_defaults

save_dir='C:\Users\mspedden\Documents\brainspine_save_bemv2';
rng(1)

if ~exist(save_dir,'dir')
    mkdir(save_dir)
end

subs = {'OP00212','OP00213',  'OP00215', 'OP00219', ...
    'OP00220', 'OP00221', 'OP00224', 'OP00225', 'OP00226'};

generic_dir = 'C:\Users\mspedden\Documents\new_leadfields_and_geom';
geomfile = fullfile(generic_dir, 'geometries_cervical_realistic.mat');

LFop='spine';
rectify=1;
fband=[10 35];
mult_comp_corr=1;

%% optional spatial smoothing (applied to observed + permutation maps before inference)
doSmooth = 1;              % 0 = no smoothing, 1 = smooth observed + permutations
smooth_fwhm_mm = 40;
smooth_radius_mm = 3*smooth_fwhm_mm;

if doSmooth
    out_suffix = sprintf('_permSmooth_%dmm', smooth_fwhm_mm);
else
    out_suffix = '';
end

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
    D=spm_eeg_load(datwithEMGmerged);
    grad_mm=D.sensors('MEG');
    ftdat = spm2fieldtrip(D);

    badchans=D.chanlabels(D.badchannels);
    cfg=[];
    cfg.channel=setdiff(ftdat.label,badchans);
    ftdat=ft_selectdata(cfg,ftdat);

    %% rectify EMG
    if rectify
        cfg=[];
        cfg.rectify='yes';
        cfg.channel='EXG1';
        ftdatr=ft_preprocessing(cfg,ftdat);
        for k=1:length(ftdat.trial)
            ftdat.trial{k}(end,:)=ftdatr.trial{k};
        end
    end

    %% load and organise spinal cord leadfields — BEM v2 with position-based channel mapping
    load('C:\Users\mspedden\Documents\bem_spine_fields\bem_v2_leadfield_cervical_realistic_bem_.mat')

    nsourcepoints = size(leadfield_cord.pos,1);
    spchanidx = find(grad_mm.coilpos(:,2) < 200);
    spchanlabs = grad_mm.label(spchanidx);

    % --- position-based channel mapping (v2 is in metres, grad_mm in mm) ---
    grad_v2 = leadfield_cord.cfg.grad;
    spc_in_v2 = zeros(numel(spchanidx), 1);
    for i = 1:numel(spchanidx)
        p = grad_mm.coilpos(spchanidx(i), :);
        d = sqrt(sum((grad_v2.coilpos*1000 - p).^2, 2));
        [~, spc_in_v2(i)] = min(d);
    end

    % clip leadfields to spinal channels in grad_mm order
    Lf = leadfield_cord;
    Lf.label = grad_mm.label(spchanidx);
    for i = 1:numel(leadfield_cord.leadfield)
        if ~isempty(leadfield_cord.leadfield{i})
            Lf.leadfield{i} = leadfield_cord.leadfield{i}(spc_in_v2, :);
        end
    end

    % dummy head model
    cfg                 = [];
    cfg.method          = 'infinite';
    cfg.siunits         = 1;
    cfg.grad            = grad_mm;
    cfg.conductivity    = 1;
    dummyvol = ft_prepare_headmodel(cfg, mesh_torso);

    %% beamforming
    cfg=[];
    cfg.output     = 'powandcsd';
    cfg.method     = 'mtmfft';
    cfg.foilim     = fband;
    cfg.tapsmofrq  = 1;
    cfg.keeptrials = 'yes';
    freqdat_tr = ft_freqanalysis(cfg, ftdat);

    cfg=[];
    cfg.avgoverfreq='yes';
    freqdat_tr = ft_selectdata(cfg, freqdat_tr);

    %% separate conditions
    statidx = find(ftdat.trialinfo==1);
    restidx  = find(ftdat.trialinfo==2);
    [nTrials,~] = min([length(statidx) length(restidx)]);

    cfg=[];
    cfg.trials = statidx(1:nTrials);
    statdat = ft_selectdata(cfg, freqdat_tr);

    cfg.trials = restidx(1:nTrials);
    restdat = ft_selectdata(cfg, freqdat_tr);

    %% mean freq data for filter construction
    cfg=[];
    cfg.output     = 'powandcsd';
    cfg.method     = 'mtmfft';
    cfg.foilim     = fband;
    cfg.tapsmofrq  = 1;
    cfg.keeptrials = 'no';
    freqdat = ft_freqanalysis(cfg, ftdat);

    cfg=[];
    cfg.avgoverfreq='yes';
    freqdat = ft_selectdata(cfg, freqdat);

    %% DICS filter
    cfg=[];
    cfg.grid                  = sources_cent;
    cfg.headmodel             = dummyvol;
    cfg.sourcemodel.leadfield = Lf;
    cfg.dics.keepfilter       = 'yes';
    cfg.dics.lambda           = 10;
    cfg.method                = 'dics';
    cfg.refchan               = 'EXG1';
    coh_source = ft_sourceanalysis(cfg, freqdat);

    %% apply filter + permutation test
    cfg=[];
    cfg.grid                  = sources_cent;
    cfg.headmodel             = dummyvol;
    cfg.sourcemodel.leadfield = Lf;
    cfg.dics.filter           = coh_source.avg.filter;
    cfg.dics.lambda           = 10;
    cfg.method                = 'dics';
    cfg.refchan               = 'EXG1';

    source_stat = ft_sourceanalysis(cfg, statdat);

    cfg.permutation    = 'yes';
    cfg.numpermutation = 500;
    source_perm = ft_sourceanalysis(cfg, statdat, restdat);

    nPerm = numel(source_perm.trialA);
    cohDiff_perm = zeros(nsourcepoints, nPerm);
    for p = 1:nPerm
        cohDiff_perm(:,p) = source_perm.trialA(p).coh - source_perm.trialB(p).coh;
    end

    coh_diff = source_perm.avgA.coh - source_perm.avgB.coh;

    %% optional spatial smoothing for inference
    if doSmooth
        Wsm = make_gaussian_smoother(sources_cent.pos, smooth_fwhm_mm, smooth_radius_mm);
        cohDiff_perm_inf = Wsm * cohDiff_perm;
        coh_diff_inf     = Wsm * coh_diff;
    else
        cohDiff_perm_inf = cohDiff_perm;
        coh_diff_inf     = coh_diff;
    end

    maxPerm = max(cohDiff_perm_inf, [], 1);
    [~, maxIdx_perm] = max(cohDiff_perm_inf, [], 1);

    %% threshold
    if mult_comp_corr
        thr95 = prctile(maxPerm, 95);
    else
        thr95 = prctile(cohDiff_perm_inf, 95, 2);
        warning('Uncorrected threshold used (per-source).')
    end

    xpos = sources_cent.pos(:,2);
    x_maxperm = xpos(maxIdx_perm);
    [~, obsMaxIdx] = max(coh_diff_inf);
    xObs = xpos(obsMaxIdx);

    %% Control plot: null maxima locations vs observed maximum
    figure('Color','w','Position',[100 100 600 450]); hold on;

    histogram(x_maxperm, 44, ...
        'FaceColor',[0.75 0.75 0.75], ...
        'EdgeColor','k', ...
        'LineWidth',0.8);

    xline(xObs, '-', ...
        'Color',[0.2 0 0], ...
        'LineWidth',2);

    xlabel('Cranio–caudal position (mm)', 'FontSize',14);
    ylabel('Count', 'FontSize',14);

    legend({'Null maxima','Observed maximum'}, ...
        'Location','best', ...
        'FontSize',14, ...
        'Box','off');

    set(gca, ...
        'FontSize',14, ...
        'LineWidth',1.2, ...
        'TickDir','out');

    box off;

    %% orientation null distribution
    W_all = coh_source.avg.filter;
    assert(size(W_all{1},1) == 3, 'Filters are not free-orientation (expected 3 x nChan).');

    A = statidx(1:nTrials);
    B = restidx(1:nTrials);
    allIdx = [A(:); B(:)];
    nA    = numel(A);
    nAll  = numel(allIdx);

    bf_labels = source_stat.cfg.channel;
    [ok, bfIdx] = ismember(bf_labels, freqdat_tr.label);
    assert(all(ok), 'Some beamformer channels not found in freqdat_tr.label');

    nPerm   = 500;
    nSource = nsourcepoints;
    ori_perm_all = nan(nSource, nPerm, 3);

    for p = 1:nPerm
        perm  = allIdx(randperm(nAll));
        permA = perm(1:nA);
        C_full = csd_from_trialset(freqdat_tr, permA);
        C_stat = real(C_full(bfIdx, bfIdx));
        for s = 1:nSource
            W = W_all{s};
            P = real(W * C_stat * W.');
            [V, D] = eig(P);
            [~, ix] = max(diag(D));
            v = V(:,ix); 
            v = v/norm(v);
            ori_perm_all(s,p,:) = v;
        end
    end
    ori_abs_all = abs(ori_perm_all);

    %% observed orientation at max significant source
    sigMask = coh_diff_inf > thr95;
    if any(sigMask)
        maskedDiff = coh_diff_inf;
        maskedDiff(~sigMask) = -inf;
        [~, obsIdx] = max(maskedDiff);
    else
        [~, obsIdx] = max(coh_diff_inf);
    end

    C_static_full = csd_from_trialset(freqdat_tr, B);
    C_static = real(C_static_full(bfIdx, bfIdx));
    Wobs = W_all{obsIdx};
    Pobs = real(Wobs * C_static * Wobs.');
    [V, D] = eig(Pobs);
    [~, ix] = max(diag(D));
    ori_obs = V(:,ix) / norm(V(:,ix));

    %% p-values + smoothed -log10(p)
    pvals = zeros(nsourcepoints,1);
    for s = 1:nsourcepoints
        permDist = cohDiff_perm_inf(s,:);
        obsVal   = coh_diff_inf(s);
        pvals(s) = (sum(permDist >= obsVal) + 1) / (nPerm + 1);
    end

    invp     = -log10(pvals);
    pthr     = 0.05;
    invpthr  = -log10(pthr);
    mask     = coh_diff_inf > thr95;
    invp_masked       = invp;
    invp_masked(~mask) = 0;

    invp_smooth = zeros(nsourcepoints,1);
    for s = 1:nsourcepoints
        permDist = sort(cohDiff_perm_inf(s,:));
        obsVal   = coh_diff_inf(s);
        xgrid    = linspace(min(permDist), max(permDist), 200);
        p_emp    = arrayfun(@(x) (sum(permDist >= x)+1)/(nPerm+1), xgrid);
        logp_smooth = smooth(xgrid, -log10(p_emp), 0.15, 'loess');
        obsVal_clamped = min(max(obsVal, xgrid(1)), xgrid(end));
        invp_smooth(s) = interp1(xgrid, logp_smooth, obsVal_clamped, 'linear');
    end

    source_p = coh_source;
    source_p.avg.coh = invp_smooth;

    cfg = [];
    cfg.parameter = 'coh';
    spine_int = ft_sourceinterpolate(cfg, source_p, mesh_wm);

    invpthr  = -log10(0.05);
    sig_vals = invp(mask & isfinite(invp));
    if isempty(sig_vals)
        clim = [invpthr invpthr+1];
    else
        clim = [invpthr max(sig_vals)];
    end

    %% clip torso mesh
    y = mesh_torso.vertices(:,2);
    keep_vert = y > -200;
    new_idx = zeros(size(keep_vert));
    new_idx(keep_vert) = 1:sum(keep_vert);
    faces_keep = all(keep_vert(mesh_torso.faces), 2);
    mesh_cut.vertices = mesh_torso.vertices(keep_vert,:);
    mesh_cut.faces    = new_idx(mesh_torso.faces(faces_keep,:));
    mesh_cut.unit     = mesh_torso.unit;

    %% plot
    ncol = 256;
    addpath('C:\Users\mspedden\Documents\fieldtrip\external\matplotlib\')
    brain_color = [0.92 0.92 0.92];
    hotmap = flipud(magma(ncol-1));
    cmap = [brain_color; hotmap];

    figure;
    cfg = [];
    cfg.figure       = 'gcf';
    cfg.method       = 'surface';
    cfg.funparameter = 'coh';
    cfg.funcolormap  = cmap;
    cfg.funcolorlim  = [2.3 2.5];
    cfg.projmethod   = 'nearest';
    cfg.surffile     = mesh_wm;
    ft_sourceplot(cfg, spine_int);
    view(-250, -1); 
    camlight;
    ax = gca; 
    ax.FontSize = 14;

    ft_plot_mesh(mesh_brain, 'facecolor', [0.8 0.3 0.3], 'facealpha', 0.07, 'edgecolor', 'none');
    ft_plot_mesh(mesh_cut,   'facecolor', [0.3 0.3 0.9], 'facealpha', 0.1,  'edgecolor', 'none'); hold on
    ft_plot_mesh(mesh_bone,  'facecolor', [0.9 0.85 0.7],'facealpha', 0.3,  'edgecolor', 'none');
    ft_plot_sens(ftdat.grad, 'coilshape', 'point', 'coilsize', 6)
    ft_plot_mesh(mesh_lungs, 'facecolor', [0.8 0.3 0.3], 'facealpha', 0.1,  'edgecolor', 'none');
    ft_plot_mesh(mesh_heart, 'facecolor', [0.8 0.3 0.3], 'facealpha', 0.1,  'edgecolor', 'none');

    %% save results
    subjResults(ss).coh_diff     = coh_diff_inf;
    subjResults(ss).thr95        = thr95;
    subjResults(ss).sig_mask     = mask;
    subjResults(ss).pos          = sources_cent.pos;
    subjResults(ss).inside       = sources_cent.inside;
    subjResults(ss).maxdiff.idx  = obsIdx;
    subjResults(ss).maxdiff.pos  = sources_cent.pos(obsIdx,:);
    subjResults(ss).maxdiff.ori  = ori_obs(:).';

end

%%-------------GROUP ANALYSIS-------------------

save(fullfile(save_dir, ['groupRes_spine_DICS_bemv2' out_suffix '.mat']), 'subjResults')

nSubjects = length(subjResults);
sig_pos   = false(nSubjects,1);

for ss = 1:nSubjects
    diffCoh = subjResults(ss).coh_diff;
    thr95   = subjResults(ss).thr95;
    if any(diffCoh > thr95)
        sig_pos(ss) = true;
    end
end

if doSmooth
    fprintf('Permutation (smoothed): %d/%d subjects show a positive effect above threshold\n', ...
        sum(sig_pos), nSubjects);
else
    fprintf('Permutation: %d/%d subjects show a positive effect above threshold\n', ...
        sum(sig_pos), nSubjects);
end

%% group prevalence
nSubs      = length(subjResults);
all_masks  = cat(2, subjResults(:).sig_mask);
group_prevalence = mean(all_masks, 2);

group_ft        = [];
group_ft.pos    = sources_cent.pos;
group_ft.inside = sources_cent.inside;
group_ft.pow    = group_prevalence;

threshold = 0.2;
group_ft.pow(group_ft.pow < threshold) = 0;

cfg = [];
cfg.parameter   = 'pow';
cfg.interpmethod = 'nearest';
group_int = ft_sourceinterpolate(cfg, group_ft, mesh_wm);

%% plot group prevalence
figure;
cfg = [];
cfg.method        = 'surface';
cfg.funparameter  = 'pow';
cfg.maskparameter = 'mask';
cfg.funcolorlim   = [threshold max(group_int.pow)];
cfg.funcolormap   = cmap;
cfg.projmethod    = 'nearest';
cfg.surffile      = mesh_wm;
cfg.opacitylim    = [threshold max(group_int.pow)];
cfg.opacitymap    = 'rampup';
ft_sourceplot(cfg, group_int);
view(-250, -1); 
camlight;
ax = gca; 
ax.FontSize = 14;

vals    = group_int.pow(:);
maxPow  = max(vals);
iMaxAll = find(vals == maxPow);
pMaxAll = group_int.pos(iMaxAll,:);

hold on
scatter3(pMaxAll(:,1), pMaxAll(:,2), pMaxAll(:,3), 80, 'k', 'filled');
ft_plot_mesh(mesh_brain, 'facecolor', [0.8 0.3 0.3], 'facealpha', 0.07, 'edgecolor', 'none');
ft_plot_mesh(mesh_cut,   'facecolor', [0.3 0.3 0.9], 'facealpha', 0.1,  'edgecolor', 'none');
ft_plot_mesh(mesh_bone,  'facecolor', [0.9 0.85 0.7],'facealpha', 0.3,  'edgecolor', 'none');
ft_plot_sens(ftdat.grad, 'coilshape', 'point', 'coilsize', 6)
ft_plot_mesh(mesh_lungs, 'facecolor', [0.8 0.3 0.3], 'facealpha', 0.1,  'edgecolor', 'none');
ft_plot_mesh(mesh_heart, 'facecolor', [0.8 0.3 0.3], 'facealpha', 0.1,  'edgecolor', 'none');

%% contiguous cluster for virtual electrode
all_masks = zeros(nsourcepoints, nSubjects);
for s = 1:nSubjects
    if isempty(subjResults(s).coh_diff)
        mask = zeros(nsourcepoints,1);
    else
        mask = subjResults(s).coh_diff > subjResults(s).thr95;
        mask = mask(:);
    end
    all_masks(:,s) = mask;
end

prevalence_loc = mean(all_masks, 2);
mask_thresh    = prevalence_loc >= threshold;
pos_thresh     = sources_cent.pos(mask_thresh,:);

distMat        = squareform(pdist(pos_thresh));
neighborRadius = 6;
adjMat         = distMat < neighborRadius;
G              = graph(adjMat);
bins           = conncomp(G);
clusterSizes   = histcounts(bins, 1:(max(bins)+1));
[~, idxMax]    = max(clusterSizes);
cluster_idx    = find(bins == idxMax);
ROIpos         = pos_thresh(cluster_idx,:);

figure; 
plot(sources_cent.pos(:,2), prevalence_loc)
hold on
for k = 1:length(ROIpos)
    plot(ROIpos(k,2), 0.2, 'r*')
end

save(fullfile(save_dir, ['cluster_spineEMG_pos_bemv2' out_suffix '.mat']), 'ROIpos')

%% visualise across subjects - 2d
nSubj = numel(subjResults);
x     = sources_cent.pos(:,2);
figure; 
hold on;

cmap = [
    27,158,119;  217,95,2;   117,112,179;
    231,41,138;  102,166,30; 230,171,2;
    166,118,29;  102,102,102; 55,126,184
    ] / 255;

for s = 1:nSubj
    cdiff = subjResults(s).coh_diff;
    thr   = subjResults(s).thr95;
    sig   = cdiff > thr;
    if any(sig), c = cmap(s,:); else, c = [0.7 0.7 0.7]; end
    for i = 1:length(x)-1
        if sig(i)&&sig(i+1)
            plot(x(i:i+1), cdiff(i:i+1), '-', 'Color', c, 'LineWidth', 1.5, 'HandleVisibility','off')
        else
            plot(x(i:i+1), cdiff(i:i+1), '-', 'Color', [0.7 0.7 0.7], 'HandleVisibility','off')
        end
    end
    plot(x(sig), cdiff(sig), '.', 'Color', c, 'MarkerSize', 12, 'HandleVisibility','off')
    h(s) = plot(nan, nan, '-', 'Color', c, 'LineWidth', 1.5);
end
yline(0, ':k', 'HandleVisibility','off')
xlabel('Cranial caudal position (mm)'); 
ylabel('Coherence difference');
if doSmooth
    title(sprintf('Significant coherence differences (smoothed inference, %d mm)', smooth_fwhm_mm))
else
    title('Significant coherence differences')
end
legend(h, arrayfun(@(s) sprintf('Participant %d', s), 1:nSubj, 'UniformOutput', false), ...
    'Location','bestoutside');
set(gca,'FontSize',13); 
grid on;

%% orientation plot
p1 = 1;
fs = 14;
nP = numel(subjResults);
sig = sig_pos(:);

oriMat = nan(nP,3);
valid  = false(nP,1);
for ss = 1:nP
    if isfield(subjResults(ss),'maxdiff') && isfield(subjResults(ss).maxdiff,'ori')
        o = subjResults(ss).maxdiff.ori(:);
        if numel(o)==3 && all(isfinite(o)) && norm(o)>0
            oriMat(ss,:) = o./norm(o);
            valid(ss)    = true;
        end
    end
end

use          = valid;
oriMat_plot  = oriMat(use,:);
idxMap       = find(use);
sig_plot     = sig(idxMap);

RL = abs(oriMat_plot(:,1));
CC = abs(oriMat_plot(:,2));
DV = abs(oriMat_plot(:,3));
X  = [RL CC DV];

figure('Color','w'); 
hold on; 
grid on;
bar(X,'grouped');
set(gca,'XTick',1:size(X,1), ...
    'XTickLabel', arrayfun(@num2str, idxMap, 'UniformOutput', false), ...
    'FontSize', fs);
xlabel('Participant','FontSize',fs);
ylabel('Axis alignment (|component|, unit vector)','FontSize',fs);
title('Orientation components at max coherence source','FontSize',fs);
legend({'Right–Left (x)','Cranial–Caudal (y)','Dorsal–Ventral (z)'}, ...
    'FontSize',fs,'Location','best');

for ii = 1:size(X,1)
    if sig_plot(ii)
        text(ii, max(X(ii,:))+0.05, '*', ...
            'HorizontalAlignment','center','FontSize',fs+4,'FontWeight','bold');
    end
end

p1_plot = find(idxMap == p1);
if ~isempty(p1_plot)
    yl = ylim;
    plot([p1_plot-0.5 p1_plot+0.5 p1_plot+0.5 p1_plot-0.5 p1_plot-0.5], ...
         [yl(1) yl(1) yl(2) yl(2) yl(1)], 'k-', 'LineWidth', 2);
end

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