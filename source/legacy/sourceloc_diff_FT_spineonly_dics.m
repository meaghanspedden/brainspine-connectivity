
% SPINE EMG COH

clear all
close all
clc

addpath('C:\Users\mspedden\Documents\brainspineconnectivity\source')
addpath('C:\Users\mspedden\Documents\spm')
spm('defaults','EEG')

addpath('C:\Users\mspedden\Documents\fieldtrip')
ft_defaults

save_dir='C:\Users\mspedden\Documents\brainspine_save';
rng(1) %for permutation testing

%n=9 for spinal cord analyses
subs = {'OP00212','OP00213',  'OP00215', 'OP00219', ...
    'OP00220', 'OP00221', 'OP00224', 'OP00225', 'OP00226'};

generic_dir = 'C:\Users\mspedden\Documents\new_leadfields_and_geom'; %where I have saved folder with brainspine leadfields and geoms (meshes)
geomfile = fullfile(generic_dir, 'geometries_cervical_realistic.mat');

LFop='spine'; %only want leadfields from spine here.
rectify=1; %EMG
fband=[10 35];
mult_comp_corr=1;

subjResults=struct();

for ss=1:length(subs)

    sub=subs{ss};

    if strcmp(sub, 'OP00224') %saved under 002
        datwithEMGmerged = fullfile('C:\Users\mspedden\Documents', ...
            ['sub-' sub], ...
            'ses-001', ...
            'meg', ...
            'pmergedoe1000mspddfflo45hi45hfcstatic_002_array1.mat');

    else

        datwithEMGmerged = fullfile('C:\Users\mspedden\Documents', ...
            ['sub-' sub], ...
            'ses-001', ...
            'meg', ...
            'pmergedoe1000mspddfflo45hi45hfcstatic_001_array1.mat');
    end

    load(geomfile)
    D=spm_eeg_load(datwithEMGmerged);
    grad_mm=D.sensors('MEG');
    ftdat = spm2fieldtrip(D);

    badchans=D.chanlabels(D.badchannels);

    %remove bad channels here.
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
            ftdat.trial{k}(end,:)=ftdatr.trial{k}; %ftdat has rectified emg
        end

    end
    %% load and organise spinal cord leadfields
    [Gx, Gy, Gz] = build_leadfield_matrices(fullfile(generic_dir,'cervical_realistic_brain_spine'), LFop);

    nsourcepoints = size(Gx,1);
    nchannels     = size(Gx,2);
    spchanidx=find(grad_mm.coilpos(:,2) < 200); %indexed locally in grad struct (same indexing as LF)
    spchanlabs=grad_mm.label(spchanidx);

    %% clip leadfields to spinal cord channels only
    Gx=Gx(:,spchanidx);
    Gy=Gy(:,spchanidx);
    Gz=Gz(:,spchanidx);

    %put leadfields into fieldtrip format
    Lf.pos    = sources_cent.pos;     % nsourcepoints x 3
    Lf.inside = sources_cent.inside;
    Lf.unit   = 'mm';
    Lf.label  = grad_mm.label(spchanidx);   % nchannels x 1 cell
    Lf.leadfielddimord = '{pos}_chan_ori';
    Lf.leadfield = cell(1,nsourcepoints);

    for k = 1:nsourcepoints
        % Combine X/Y/Z components like FT is used to
        Lf.leadfield{k} = [Gx(k,:)' Gy(k,:)' Gz(k,:)']; % nchannels x 3
    end

    % 2. dummy head model for input config only (not actually used)
    cfg                     = [];
    cfg.method              = 'infinite';
    cfg.siunits=1;
    cfg.grad=grad_mm;
    cfg.conductivity = 1;

    dummyvol = ft_prepare_headmodel(cfg,mesh_torso);


    %% beamforming----------------------------
    %1. get trial wise freq dat
    cfg=[];
    cfg.output     = 'powandcsd';
    cfg.method     = 'mtmfft';
    cfg.foilim     = fband;
    cfg.tapsmofrq  = 1;
    cfg.keeptrials='yes';
    freqdat_tr=ft_freqanalysis(cfg,ftdat);

    cfg=[];
    cfg.avgoverfreq='yes';
    freqdat_tr=ft_selectdata(cfg,freqdat_tr);%need to average over frequency - required input to permutation test

    %% separate conditions

    statidx=find(ftdat.trialinfo==1);
    restidx=find(ftdat.trialinfo==2);

    [nTrials,~] = min([length(statidx) length(restidx)]); %ensure same n trials across conditions

    cfg=[];
    cfg.trials=statidx(1:nTrials);
    statdat=ft_selectdata(cfg,freqdat_tr); %trial wise data pr condition

    cfg.trials=restidx(1:nTrials);
    restdat=ft_selectdata(cfg,freqdat_tr);

    %% get struct with mean freq data for constructing bf filter
    cfg=[];
    cfg.output     = 'powandcsd';
    cfg.method     = 'mtmfft';
    cfg.foilim     = fband;
    cfg.tapsmofrq  = 1;
    cfg.keeptrials = 'no';
    freqdat=ft_freqanalysis(cfg,ftdat);

    cfg=[];
    cfg.avgoverfreq='yes';
    freqdat=ft_selectdata(cfg,freqdat);

    %% DICS filter constructed on both conditions

    cfg=[];
    cfg.grid = sources_cent;
    cfg.headmodel=dummyvol;
    cfg.sourcemodel.leadfield=Lf;
    cfg.dics.keepfilter='yes';
    cfg.dics.lambda=10;
    cfg.method = 'dics';
    cfg.refchan='EXG1'; %EMG

    coh_source = ft_sourceanalysis(cfg,freqdat);

    %% apply to each condition separately, do permutation test
    cfg=[];
    cfg.grid = sources_cent;
    cfg.headmodel=dummyvol;
    cfg.sourcemodel.leadfield=Lf;
    cfg.dics.filter=coh_source.avg.filter;
    cfg.dics.lambda=10;
    cfg.method = 'dics';
    cfg.refchan='EXG1';

    source_stat = ft_sourceanalysis(cfg,statdat);
    %source_rest = ft_sourceanalysis(cfg,restdat);

    cfg.permutation = 'yes';
    cfg.numpermutation=500;
    source_perm = ft_sourceanalysis(cfg, statdat, restdat); %permutation test

    nPerm = numel(source_perm.trialA);
    cohDiff_perm = zeros(nsourcepoints, nPerm);

    for p = 1:nPerm
        cohDiff_perm(:,p) = source_perm.trialA(p).coh - source_perm.trialB(p).coh; % A-B
    end

    % observed (unpermuted) coherence difference
    coh_diff = source_perm.avgA.coh - source_perm.avgB.coh;  % nSourcePoints x 1

    % max over sources per permutation (for global threshold)
    maxPerm = max(cohDiff_perm, [], 1);              % 1 x nPerm
    [~, maxIdx_perm] = max(cohDiff_perm, [], 1);     % 1 x nPerm (argmax location)

    %% 2) Global threshold (or uncorrected per-source threshold)

    if mult_comp_corr
        thr95 = prctile(maxPerm, 95);        % scalar global threshold
    else
        thr95 = prctile(cohDiff_perm, 95, 2); % nSourcePoints x 1 (uncorrected)
        warning('Uncorrected threshold used (per-source).')
    end

   %% 3) Control plot: where null maxima land + probability they fall in ROI

nsourcepoints = numel(xpos);

roi_idx = 16:25;      % expected ROI (C6–T1)
sig_idx = 16:18;      % observed significant run for sub-annotation (optional)

% --- Null max locations (global max index per perm) ---
x_maxperm = xpos(maxIdx_perm);

% --- Observed global max location ---
[~, obsMaxIdx] = max(coh_diff);
xObs = xpos(obsMaxIdx);

% --- ROI bounds in mm ---
x_roi_lo = xpos(min(roi_idx));
x_roi_hi = xpos(max(roi_idx));

x_sig_lo = xpos(min(sig_idx));
x_sig_hi = xpos(max(sig_idx));

% --- Location (descriptive) probabilities ---
p_loc_roi = (sum(ismember(maxIdx_perm, roi_idx)) + 1) / (numel(maxIdx_perm) + 1);
p_loc_sig = (sum(ismember(maxIdx_perm, sig_idx)) + 1) / (numel(maxIdx_perm) + 1);

% --- ROI max (inferential) statistic ---
obs_roi_max  = max(coh_diff(roi_idx));
null_roi_max = max(cohDiff_perm(roi_idx,:), [], 1);
p_roi = (sum(null_roi_max >= obs_roi_max) + 1) / (nPerm + 1);

%% Plot
figure('Color','w','Position',[100 100 650 450]); hold on;

h = histogram(x_maxperm, 44, ...
    'FaceColor',[0.75 0.75 0.75], ...
    'EdgeColor','k', ...
    'LineWidth',0.8);

yl = ylim;

% Shade expected ROI (light)
pROI = patch([x_roi_lo x_roi_hi x_roi_hi x_roi_lo], ...
             [yl(1)    yl(1)    yl(2)    yl(2)], ...
             [0 0 0], 'FaceAlpha',0.08, 'EdgeColor','none');

% Shade observed sig band (darker) - optional
pSIG = patch([x_sig_lo x_sig_hi x_sig_hi x_sig_lo], ...
             [yl(1)    yl(1)    yl(2)    yl(2)], ...
             [0 0 0], 'FaceAlpha',0.16, 'EdgeColor','none');

% Keep bars on top of shading
uistack(h,'top');

% Observed maximum line
%lObs = xline(xObs, '-', 'Color',[0.2 0 0], 'LineWidth',2);

% ROI boundary lines (optional)
xline(x_roi_lo, '--', 'Color',[0 0 0], 'LineWidth',1.2);
xline(x_roi_hi, '--', 'Color',[0 0 0], 'LineWidth',1.2);

% Text block (make ROI-max the headline)
txt = sprintf([ ...
    'ROI-max test within C6–T1 (16–25):  p = %.4f\n' ...
    'Peak-location (descriptive):  P(null peak in ROI) = %.3f\n' ],...
    p_roi, p_loc_roi);

% txt = sprintf([ ...
%     'ROI-max test within C6–T1 (16–25):  p = %.4f\n' ...
%     'Peak-location (descriptive):  P(null peak in ROI) = %.3f\n' ...
%     'Peak-location (descriptive):  P(null peak in 16–18) = %.3f' ], ...
%     p_roi, p_loc_roi, p_loc_sig);

% text(min(xlim) + 0.02*range(xlim), 115, txt, ...
%     'HorizontalAlignment','left', ...
%     'VerticalAlignment','top', ...
%     'FontSize',12);

xlabel('Cranio–caudal position (mm)', 'FontSize',14);
ylabel('Count', 'FontSize',14);

legend([h pROI pSIG], ...
    {'Null peak locations','Expected ROI (C6–T1)','Observed significant'}, ...
    'Location','best', 'FontSize',12, 'Box','off');

set(gca, 'FontSize',14, 'LineWidth',1.2, 'TickDir','out');
box off;


%% Orientation null at EACH sourcepoint (not at peak)

W_all = coh_source.avg.filter;   % cell{nSource}, each is 3 x nChan
assert(size(W_all{1},1) == 3, 'Filters are not free-orientation (expected 3 x nChan).');

% Balanced trial pools (as you had)
A = statidx(1:nTrials);
B = restidx(1:nTrials);
allIdx = [A(:); B(:)];
nA = numel(A);
nAll = numel(allIdx);

% Beamformer channel labels ordering
bf_labels = source_stat.cfg.channel;
[ok, bfIdx] = ismember(bf_labels, freqdat_tr.label);
assert(all(ok), 'Some beamformer channels not found in freqdat_tr.label');

nPerm = 500;  % or cfg.numpermutation / size you want
nSource = nsourcepoints;

ori_perm_all = nan(nSource, nPerm, 3);   % source x perm x component

for p = 1:nPerm

    % permute trial labels and take "static" pool (as before)
    perm  = allIdx(randperm(nAll));
    permA = perm(1:nA);

    % compute CSD once per permutation
    C_full = csd_from_trialset(freqdat_tr, permA);
    C_stat = real(C_full(bfIdx, bfIdx));

    % now compute orientation at ALL sources for this permutation
    for s = 1:nSource
        W = W_all{s};                 % 3 x nChan
        P = real(W * C_stat * W.');   % 3 x 3

        [V, D] = eig(P);
        [~, ix] = max(diag(D));
        v = V(:,ix);
        v = v / norm(v);

        ori_perm_all(s,p,:) = v;
    end
end

% Because of sign ambiguity, usually work with absolute components:
ori_abs_all = abs(ori_perm_all);   % nSource x nPerm x 3

    %% 5) Observed (real) orientation at max SIGNIFICANT coherence-difference location

    sigMask = coh_diff > thr95;                 % scalar thr95 or per-source vector thr95

    if any(sigMask)
        maskedDiff = coh_diff;
        maskedDiff(~sigMask) = -inf;
        [~, obsIdx] = max(maskedDiff);          % max significant diff
    else
        [~, obsIdx] = max(coh_diff);            % fallback
    end

    % Static/rest CSD from true rest trials (B)
    C_static_full = csd_from_trialset(freqdat_tr, B);
    C_static = real(C_static_full(bfIdx, bfIdx));

    Wobs = W_all{obsIdx};
    Pobs = real(Wobs * C_static * Wobs.');

    [V, D] = eig(Pobs);
    [~, ix] = max(diag(D));
    ori_obs = V(:,ix) / norm(V(:,ix));          % 3 x 1

%     if ~any(sigMask)
%         ori_obs = [NaN NaN NaN];
%     end

    % resolve sign ambiguity for plotting/comparison
    for p = 1:nPerm
        if dot(ori_perm(p,:), ori_obs) < 0
            ori_perm(p,:) = -ori_perm(p,:);
        end
    end

    %% 6) Summary plot: axis components (null vs observed)

figure('Color','w','Position',[750 100 600 450]); hold on;

% Mean absolute permuted orientation components
b = bar(mean(abs(ori_perm), 1), ...
    'FaceColor',[0.75 0.75 0.75], ...
    'EdgeColor','k', ...
    'LineWidth',0.8);

% Observed orientation components
plot(1:3, abs(ori_obs), 'o', ...
    'MarkerSize',10, ...
    'MarkerEdgeColor',[0.2 0 0], ...
    'MarkerFaceColor',[0.2 0 0], ...
    'LineWidth',1.5);

% Axes
set(gca, ...
    'XTick',1:3, ...
    'XTickLabel',{'L–R','C–C','D–V'}, ...
    'FontSize',14, ...
    'LineWidth',1.2, ...
    'TickDir','out');

ylabel('|Orientation component|', 'FontSize',14);

legend({'Permuted (mean)','Observed'}, ...
    'Location','best', ...
    'FontSize',14, ...
    'Box','off');

box off;

    %% 7) One-sided permutation p-values + -log10(p) masked map (for plotting)

    pvals = zeros(nsourcepoints, 1);
    for s = 1:nsourcepoints
        permDist = cohDiff_perm(s, :);
        obsVal   = coh_diff(s);
        pvals(s) = (sum(permDist >= obsVal) + 1) / (nPerm + 1);
    end

    invp = -log10(pvals);
    pthr = 0.05;
    invpthr = -log10(pthr);

    % Apply same significance mask as thresholding
    mask = coh_diff > thr95;
    invp_masked = invp;
    invp_masked(~mask) = 0;  % or NaN

    % Put into a source structure for interpolation/plotting
    source_p = coh_source;
    source_p.avg.coh = invp_masked;
    
   source_p.avg.coh =  invp_smooth;

    cfg = [];
    cfg.parameter = 'coh';
    spine_int = ft_sourceinterpolate(cfg, source_p, mesh_wm);
    % --- robust color limits ---
invpthr = -log10(0.05);  % or whatever you use
sig_vals = invp(mask & isfinite(invp));

if isempty(sig_vals)
    % Nothing significant: choose a harmless range so ft_sourceplot doesn't crash
    clim = [invpthr invpthr+1];  % any non-zero width range is fine
else
    clim = [invpthr max(sig_vals)];
end

    %% 8) Clip torso mesh (for cleaner plotting)

    y = mesh_torso.vertices(:,2);
    keep_vert = y > -200;

    new_idx = zeros(size(keep_vert));
    new_idx(keep_vert) = 1:sum(keep_vert);

    faces_keep = all(keep_vert(mesh_torso.faces), 2);
    mesh_cut.vertices = mesh_torso.vertices(keep_vert,:);
    mesh_cut.faces    = new_idx(mesh_torso.faces(faces_keep,:));
    mesh_cut.unit     = mesh_torso.unit;

    %% 9) Plot -log10(p) on spinal mesh

    ncol = 256;
    addpath('C:\Users\mspedden\Documents\fieldtrip\external\matplotlib\')
    brain_color = [0.92 0.92 0.92];
    hotmap = flipud(magma(ncol-1));

    cmap = [brain_color; hotmap];
 
    figure;  
    cfg = [];
    cfg.figure      = 'gcf';
    cfg.method      = 'surface';
    cfg.funparameter= 'coh';
    cfg.funcolormap = cmap;
    cfg.funcolorlim = [2.3 2.5];
    cfg.projmethod  = 'nearest';
    cfg.surffile    = mesh_wm;
    ft_sourceplot(cfg, spine_int);

    view(-250, -1);
    camlight;
    ax = gca;
    ax.FontSize = 14;
     ft_plot_mesh(mesh_brain, 'facecolor', [0.8 0.3 0.3], 'facealpha', 0.07, 'edgecolor', 'none');
    ft_plot_mesh(mesh_cut, 'facecolor', [0.3 0.3 0.9], 'facealpha', 0.1, 'edgecolor', 'none'); hold on
    ft_plot_mesh(mesh_bone, 'facecolor', [0.9 0.85 0.7], 'facealpha', 0.3, 'edgecolor', 'none');
    ft_plot_sens(ftdat.grad,'coilshape','point','coilsize',6)
    ft_plot_mesh(mesh_lungs, 'facecolor', [0.8 0.3 0.3], 'facealpha', 0.1, 'edgecolor', 'none');
    ft_plot_mesh(mesh_heart, 'facecolor', [0.8 0.3 0.3], 'facealpha', 0.1, 'edgecolor', 'none');


    %% save results

subjResults(ss).coh_diff = coh_diff;          % A - B coherence difference (nSourcePoints x 1)
subjResults(ss).thr95    = thr95;             % significance threshold

% significance / geometry 
subjResults(ss).sig_mask = mask;              % logical mask in source space
subjResults(ss).pos      = sources_cent.pos;  % source positions
subjResults(ss).inside   = sources_cent.inside;

% max significant source
subjResults(ss).maxdiff.idx = obsIdx;                         % source index
subjResults(ss).maxdiff.pos = sources_cent.pos(obsIdx,:);     % OPTIONAL but useful
subjResults(ss).maxdiff.ori = ori_obs(:).';                   % 1 x 3 orientation
end

%%-------------GROUP ANALYSIS-------------------
%%----------------------------------------------


save('groupRes_spine_DICS.mat', 'subjResults')

nSubjects = length(subjResults);
sig_pos = false(nSubjects,1);

%number of subjects with at least one significant source anywhere
for ss = 1:nSubjects
    diffCoh = subjResults(ss).coh_diff;   % source_perm.avgA.coh - avgB.coh
    thr95    = subjResults(ss).thr95; % 95th percentile from permutation

    if any(diffCoh > thr95)
        sig_pos(ss) = true;
    end
end

fprintf('Permutation: %d/%d subjects show a positive effect above threshold\n', ...
    sum(sig_pos), nSubjects);


%% group prevalence

nSubs = length(subjResults);

all_masks = cat(2, subjResults(:).sig_mask);

% Compute prevalence
group_prevalence = mean(all_masks, 2);

% Create a group source structure for plotting
group_source = subjResults(1); % use one subject as template
group_source.pow = group_prevalence;
group_source = rmfield(group_source, {'coh_diff', 'sig_mask'}); % clean up

group_ft = [];
group_ft.pos = group_source.pos;
group_ft.inside = group_source.inside;
group_ft.pow = group_prevalence;

%% Interpolate group map onto the mesh
threshold = 0.2; 
group_ft.pow(group_ft.pow < threshold) = 0;  % threshold source points

cfg = [];
cfg.parameter = 'pow';
cfg.interpmethod = 'nearest';
group_int = ft_sourceinterpolate(cfg, group_ft, mesh_wm);


%% Plot group prevalence map
figure;

cfg = [];
cfg.method = 'surface';
cfg.funparameter = 'pow';
cfg.maskparameter = 'mask';
cfg.funcolorlim =  [threshold max(group_int.pow)];
cfg.funcolormap = cmap;
cfg.projmethod = 'nearest';
cfg.surffile = mesh_wm;
cfg.opacitylim    = [threshold max(group_int.pow)];
cfg.opacitymap    = 'rampup';
ft_sourceplot(cfg, group_int);
view(-250, -1);
camlight;
ax = gca;
ax.FontSize = 14;
hold on
ft_plot_mesh(mesh_brain, 'facecolor', [0.8 0.3 0.3], 'facealpha', 0.07, 'edgecolor', 'none');
ft_plot_mesh(mesh_cut, 'facecolor', [0.3 0.3 0.9], 'facealpha', 0.1, 'edgecolor', 'none'); hold on
ft_plot_mesh(mesh_bone, 'facecolor', [0.9 0.85 0.7], 'facealpha', 0.3, 'edgecolor', 'none');
ft_plot_sens(ftdat.grad,'coilshape','point','coilsize',6)
ft_plot_mesh(mesh_lungs, 'facecolor', [0.8 0.3 0.3], 'facealpha', 0.1, 'edgecolor', 'none');
ft_plot_mesh(mesh_heart, 'facecolor', [0.8 0.3 0.3], 'facealpha', 0.1, 'edgecolor', 'none');
hold on

%% Find contiguous cluster for virtual electrode
all_masks = zeros(nsourcepoints, nSubjects);

for s = 1:nSubjects
    if isempty(subjResults(s).coh_diff)
        mask = zeros(nsourcepoints,1);   % no significant sources
    else
        mask = subjResults(s).coh_diff > subjResults(s).thr95;
        mask = mask(:);
    end
    all_masks(:,s) = mask;
end

prevalence_loc = mean(all_masks, 2);  % nSourcePoints x 1

mask_thresh = prevalence_loc >= threshold;

% Extract positions of thresholded sources
pos_thresh = sources_cent.pos(mask_thresh,:);

% Compute pairwise distances
distMat = squareform(pdist(pos_thresh));

% Define adjacency points within 6 mm as neighbours
neighborRadius = 6;
adjMat = distMat < neighborRadius;

G = graph(adjMat);
bins = conncomp(G);

% Compute cluster sizes
clusterSizes = histcounts(bins, 1:(max(bins)+1));

% Get largest contiguous cluster
[~, idxMax] = max(clusterSizes);

% Indices of points in largest cluster (in pos_thresh)
cluster_idx = find(bins == idxMax);

% Positions of largest cluster in source space
ROIpos = pos_thresh(cluster_idx,:);

%Plot positions of largest cluster
figure; plot(sources_cent.pos(:,2), prevalence_loc)
hold on
for k=1:length(ROIpos)
    plot(ROIpos(k,2), 0.2, 'r*')
end

save(fullfile(save_dir, 'cluster_spineEMG_pos.mat'), 'ROIpos')

%% binomial p and CI
x = sum(sig_pos);
n = nSubjects;
alpha=0.05;
phat = x/n;
lower = betainv(alpha/2, x,     n-x+1);
upper = betainv(1-alpha/2, x+1, n-x);

ci = [lower upper];
p = 1 - binocdf(x-1, n, alpha);


% figure;
% scatter3(group_ft.pos(:,1), group_ft.pos(:,2), group_ft.pos(:,3), ...
%     50, group_prevalence, 'filled');
% colorbar;
% caxis([0 1]);  % prevalence between 0 and 1
% axis equal;
% xlabel('X'); ylabel('Y'); zlabel('Z');
% title('Group Prevalence (no interpolation)');



%% visualise across subjects - 2d

nSubj=numel(subjResults); x=sources_cent.pos(:,2); figure; hold on;

cmap = [
    27,158,119;
    217,95,2;
    117,112,179;
    231,41,138;
    102,166,30;
    230,171,2;
    166,118,29;
    102,102,102;
    55,126,184
    ] / 255;

for s=1:nSubj
    cdiff=subjResults(s).coh_diff; thr=subjResults(s).thr95; sig=cdiff>thr;
    if any(sig), c=cmap(s,:); else, c=[0.7 0.7 0.7]; end
    for i=1:length(x)-1
        if sig(i)&&sig(i+1)
            plot(x(i:i+1),cdiff(i:i+1),'-','Color',c,'LineWidth',1.5,'HandleVisibility','off')
        else
            plot(x(i:i+1),cdiff(i:i+1),'-','Color',[0.7 0.7 0.7],'HandleVisibility','off')
        end
    end
    plot(x(sig),cdiff(sig),'.','Color',c,'MarkerSize',12,'HandleVisibility','off')
    h(s)=plot(nan,nan,'-','Color',c,'LineWidth',1.5);
end
yline(0,':k','HandleVisibility','off')
xlabel('Cranial caudal position (mm)'); ylabel('Coherence difference'); title('Significant coherence differences')
legend(h, arrayfun(@(s) sprintf('Participant %d', s), 1:nSubj, 'UniformOutput', false), ...
    'Location', 'bestoutside');

set(gca, 'FontSize', 13)
grid on;

%% plot orientation max coherence across subjects 
p1 = 1;
fs = 14;

nP = numel(subjResults);

% --- significance mask ---
sig = sig_pos(:);   % logical nP x 1

% --- collect orientations ---
oriMat = nan(nP,3);
valid  = false(nP,1);

for ss = 1:nP
    if isfield(subjResults(ss),'maxdiff') && isfield(subjResults(ss).maxdiff,'ori')
        o = subjResults(ss).maxdiff.ori(:);
        if numel(o)==3 && all(isfinite(o)) && norm(o)>0
            oriMat(ss,:) = o./norm(o);
            valid(ss) = true;
        end
    end
end

% --- keep all VALID participants ---
use = valid;
oriMat_plot = oriMat(use,:);
idxMap = find(use);              % plotted index -> original participant index
sig_plot = sig(idxMap);          % significance for plotted subjects

% components
RL = abs(oriMat_plot(:,1));
CC = abs(oriMat_plot(:,2));
DV = abs(oriMat_plot(:,3));
X = [RL CC DV];

% --- plot ---
figure('Color','w'); hold on; grid on;
barHandles = bar(X,'grouped');

nPlot = size(X,1);

set(gca,'XTick',1:nPlot, ...
        'XTickLabel', arrayfun(@num2str, idxMap, 'UniformOutput', false), ...
        'FontSize', fs);

xlabel('Participant','FontSize',fs);
ylabel('Axis alignment (|component|, unit vector)','FontSize',fs);
title('Orientation components at max coherence source','FontSize',fs);

legend({'Right–Left (x)','Cranial–Caudal (y)','Dorsal–Ventral (z)'}, ...
       'FontSize',fs, 'Location','best');

% --- add significance stars ---
yl = ylim;
starY = yl(2) * 0.95;   % position near top

for ii = 1:nPlot
    if sig_plot(ii)
        yMax = max(X(ii,:));
        text(ii, yMax + 0.05, '*', ...
            'HorizontalAlignment','center', ...
            'FontSize',fs+4, ...
            'FontWeight','bold');
    end
end

% --- highlight Participant 1 if present ---
p1_plot = find(idxMap == p1);
if ~isempty(p1_plot)
    yl = ylim;
    plot([p1_plot-0.5 p1_plot+0.5 p1_plot+0.5 p1_plot-0.5 p1_plot-0.5], ...
         [yl(1) yl(1) yl(2) yl(2) yl(1)], 'k-', 'LineWidth', 2);
end


out_spine = plot_bayesprev_posterior(sig_pos, 0.05);
%% sorted by height

% heighttable=readtable('C:\Users\mspedden\Documents\SC_subs_heights.csv');
% heights=heighttable.Var2;
%
% [sortedHeights, sortIdx] = sort(heights, 'descend');  % tallest first
% subjResultsSorted = subjResults(sortIdx);
%
% cmapSorted = cmap(sortIdx, :);
%
% nSubj = numel(subjResultsSorted);
% figure; hold on;
%
% for s = 1:nSubj
%     cdiff = subjResultsSorted(s).coh_diff;
%     thr   = subjResultsSorted(s).thr95;
%     sig   = cdiff > thr;
%
%     if any(sig)
%         c = cmapSorted(s,:);
%     else
%         c = [0.7 0.7 0.7];
%     end
%
%     for i = 1:length(x)-1
%         if sig(i) && sig(i+1)
%             plot(x(i:i+1), cdiff(i:i+1), '-', 'Color', c, 'LineWidth', 1.5, 'HandleVisibility', 'off')
%         else
%             plot(x(i:i+1), cdiff(i:i+1), '-', 'Color', [0.7 0.7 0.7], 'HandleVisibility', 'off')
%         end
%     end
%
%     plot(x(sig), cdiff(sig), '.', 'Color', c, 'MarkerSize', 12, 'HandleVisibility', 'off')
%     h(s) = plot(nan, nan, '-', 'Color', c, 'LineWidth', 1.5);
% end
%
% yline(0, ':k', 'HandleVisibility', 'off');
% xlabel('Cranial caudal position (mm)');
% ylabel('Coherence difference');
% title('Significant coherence differences height sorted');
% legend(h, arrayfun(@(s) sprintf('Subj %d', s), 1:nSubj, 'UniformOutput', false), 'Location', 'bestoutside');
% set(gca, 'FontSize', 13)
% grid on;



