%% step_B_M1VE_spine_coherence.m
% DICS coherence between M1 virtual electrode and spinal cord sources.
% Contraction vs rest permutation test. All subjects in both lists.

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
lf_path        = 'C:\Leadfields meshes\leadfield_experimental_bslaw_experimental.mat';

subs_brain = {'OP00212','OP00213','OP00215','OP00219', ...
              'OP00225','OP00221','OP00224'};
subs_spine = {'OP00212','OP00213','OP00215','OP00219', ...
              'OP00220','OP00221','OP00224','OP00225','OP00226'};
subs_both  = intersect(subs_brain, subs_spine, 'stable');

fband          = [10 35];
numpermutation = 500;
mult_comp_corr = 1;
lambda         = '10%';
fwhm_mm        = 20;
radius_mm      = 3 * (fwhm_mm / 2.355);
brain_suffix   = '_brain_pct10';
out_suffix     = '_M1VE_spine_true_bem';
roi_idx        = 14:27;   % C6-T1 anatomical ROI for plots
rng(1);

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

ncol     = 256;
addpath(fullfile(fieldtrip_path,'external','matplotlib'));
cmap_hot = [[0.92 0.92 0.92]; flipud(magma(ncol-1))];

%% =========================================================================
%  LOAD SPINE GEOMETRY
%% =========================================================================
fprintf('Loading geometry...\n');
geom_exp     = load(geomfile_spine);
sources_cent = geom_exp.sources_cent;
mesh_torso   = geom_exp.mesh_torso;
mesh_wm      = geom_exp.mesh_wm;
mesh_wm.unit = 'mm';

cord_pos      = sources_cent.pos(:,2);
nsourcepoints = size(sources_cent.pos, 1);

%% =========================================================================
%  BUILD SMOOTHER ONCE
%% =========================================================================
fprintf('Building Gaussian smoother (FWHM=%d mm)...\n', fwhm_mm);
Wsm = make_gaussian_smoother(sources_cent.pos, fwhm_mm, radius_mm);

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

    %% Load brain VE
    ve_file = fullfile(save_dir, sprintf('sub%s_VE_brain_M1%s.mat', sub, brain_suffix));
    assert(exist(ve_file,'file')==2, 'Brain VE not found: %s', ve_file);
    loaded_ve  = load(ve_file, 'VE_brain');
    VE_brain   = loaded_ve.VE_brain;
    fprintf('  Brain VE loaded\n');

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

    %% Build spine leadfield
    fprintf('  Building spine leadfield...\n');
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

    %% Append brain VE to spine data
    spinedat            = ftdat;
    trialinfo           = ftdat.trialinfo;
    spinedat            = ft_appenddata([], spinedat, VE_brain);
    spinedat.trialinfo  = trialinfo;

    %% Volume conductor
    cfg = []; cfg.method = 'infinite'; cfg.siunits = 1;
    cfg.grad = grad_mm; cfg.conductivity = 1;
    dummyvol = ft_prepare_headmodel(cfg, mesh_torso);

    %% Frequency data
    cfg_sel = []; cfg_sel.channel = [Lf.label; {'virtualchannel001'}];
    cfg_fr  = []; cfg_fr.output = 'powandcsd'; cfg_fr.method = 'mtmfft';
    cfg_fr.foilim = fband; cfg_fr.tapsmofrq = 1; cfg_fr.keeptrials = 'yes';
    cfg_av  = []; cfg_av.avgoverfreq = 'yes';

    freqdat_tr = ft_freqanalysis(cfg_fr, spinedat);
    freqdat_tr = ft_selectdata(cfg_av,  freqdat_tr);

    statidx = find(trialinfo == 1);
    restidx = find(trialinfo == 2);
    nTrials = min(numel(statidx), numel(restidx));

    cfg = []; cfg.trials = statidx(1:nTrials);
    statdat = ft_selectdata(cfg, freqdat_tr);
    cfg = []; cfg.trials = restidx(1:nTrials);
    restdat = ft_selectdata(cfg, freqdat_tr);

    % Combined for common spatial filter
    cfg_fr2 = []; cfg_fr2.output = 'powandcsd'; cfg_fr2.method = 'mtmfft';
    cfg_fr2.foilim = fband; cfg_fr2.tapsmofrq = 1; cfg_fr2.keeptrials = 'no';
    freqdat = ft_freqanalysis(cfg_fr2, spinedat);
    cfg_av2 = []; cfg_av2.avgoverfreq = 'yes';
    freqdat = ft_selectdata(cfg_av2, freqdat);

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

    %% Common spatial filter
    fprintf('  Computing common spatial filter...\n');
    cfg_dics = [];
    cfg_dics.sourcemodel     = sourcemodel;
    cfg_dics.headmodel       = dummyvol;
    cfg_dics.dics.keepfilter = 'yes';
    cfg_dics.dics.lambda     = lambda;
    cfg_dics.method          = 'dics';
    cfg_dics.refchan         = 'virtualchannel001';
    coh_source = ft_sourceanalysis(cfg_dics, freqdat);

    %% Permutation test
    fprintf('  Running permutation test (%d permutations)...\n', numpermutation);
    cfg_perm = [];
    cfg_perm.sourcemodel    = sourcemodel;
    cfg_perm.headmodel      = dummyvol;
    cfg_perm.dics.filter    = coh_source.avg.filter;
    cfg_perm.dics.lambda    = lambda;
    cfg_perm.method         = 'dics';
    cfg_perm.refchan        = 'virtualchannel001';
    cfg_perm.permutation    = 'yes';
    cfg_perm.numpermutation = numpermutation;
    source_perm = ft_sourceanalysis(cfg_perm, statdat, restdat);

    nPerm = numel(source_perm.trialA);
    [coh_diff, cohDiff_perm] = extract_coh_diff(source_perm, nsourcepoints, nPerm);

    %% Smoothing
    cohDiff_perm = Wsm * cohDiff_perm;
    coh_diff     = Wsm * coh_diff;

    %% Threshold
    thr95 = compute_threshold(cohDiff_perm, mult_comp_corr, nsourcepoints);
    mask  = coh_diff > thr95;

    pvals       = (sum(cohDiff_perm >= coh_diff, 2) + 1) / (numpermutation + 1);
    invp        = -log10(pvals);
    invp_smooth = smooth_invp(coh_diff, cohDiff_perm, nsourcepoints, nPerm);
    invpthr     = -log10(0.05);

    fprintf('  Threshold: %.6f  |  Sig sources: %d / %d\n', ...
        thr95, sum(mask), nsourcepoints);
    [peak_val, peak_idx] = max(coh_diff);
    fprintf('  Peak coh_diff: %.4e at y=%.1f mm (source %d)\n', ...
        peak_val, cord_pos(peak_idx), peak_idx);
    fprintf('  Subject time: %.1f min\n', toc(t_sub)/60);

    %% Interpolate onto spine mesh
    source_p         = coh_source;
    source_p.avg.coh = invp_smooth;
    cfg_interp = []; cfg_interp.parameter = 'coh';
    spine_int  = ft_sourceinterpolate(cfg_interp, source_p, mesh_wm);

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
    cfg_plot.funcolorlim  = clim_spine; cfg_plot.projmethod = 'nearest';
    cfg_plot.surffile = mesh_wm;
    ft_sourceplot(cfg_plot, spine_int);
    view(-250,-1); camlight; ax = gca; ax.FontSize = 14;
    hold on;
    ft_plot_mesh(mesh_wm,'facecolor',[0.7 0.7 0.7],'facealpha',0.1,'edgecolor','none');
    ft_plot_mesh(mesh_cut,'facecolor',[0.3 0.3 0.9],'facealpha',0.1,'edgecolor','none');
    title(sprintf('%s — M1VE-spine coherence (true BEM)', sub),'Interpreter','none');
    savefig(hfig_spine, fullfile(fig_dir, ...
        sprintf('stepB_sub%s_M1VE_spine%s.fig', sub, out_suffix)));
    close(hfig_spine);

    %% Null maxima histogram
    [~, maxIdx_perm] = max(cohDiff_perm, [], 1);
    [~, obsMaxIdx]   = max(coh_diff);

    hfig_null = figure('Color','w','Position',[100 100 600 450]);
    hold on;
    histogram(cord_pos(maxIdx_perm), 44, ...
        'FaceColor',[0.75 0.75 0.75],'EdgeColor','k','LineWidth',0.8);
    xline(cord_pos(obsMaxIdx),'-','Color',[0.2 0 0],'LineWidth',2);
    xlabel('Cranio-caudal position (mm)','FontSize',14);
    ylabel('Count','FontSize',14);
    legend({'Null maxima','Observed maximum'},'Location','best','FontSize',14,'Box','off');
    set(gca,'FontSize',14,'LineWidth',1.2,'TickDir','out'); box off;
    title(sprintf('%s — null maxima (M1VE-spine)', sub),'Interpreter','none');
    savefig(hfig_null, fullfile(fig_dir, ...
        sprintf('stepB_sub%s_null_maxima%s.fig', sub, out_suffix)));
    close(hfig_null);

    %% Store
    subjResults(ss).subjID      = sub;
    subjResults(ss).coh_diff    = coh_diff;
    subjResults(ss).thr95       = thr95;
    subjResults(ss).sig_mask    = mask;
    subjResults(ss).pos         = sources_cent.pos;
    subjResults(ss).inside      = logical(Lf.inside);
    subjResults(ss).invp_smooth = invp_smooth;
    subjResults(ss).cord_pos    = cord_pos;

    save(fullfile(save_dir, sprintf('stepB_result_sub%s%s.mat', sub, out_suffix)), ...
        'coh_diff','cohDiff_perm','thr95','mask','invp_smooth','cord_pos');
end

%% Save group
save(fullfile(save_dir, ['groupRes_M1VE_spine' out_suffix '.mat']), 'subjResults');
fprintf('\nAll subjects complete. Running group analysis...\n');

%% =========================================================================
%  GROUP ANALYSIS
%% =========================================================================
nSubjects     = numel(subjResults);
all_masks     = zeros(nsourcepoints, nSubjects);
for s = 1:nSubjects
    all_masks(:,s) = double(subjResults(s).coh_diff > subjResults(s).thr95);
end

sig_pos = false(nSubjects,1);
for s = 1:nSubjects
    sig_pos(s) = any(all_masks(:,s));
end
fprintf('  %d/%d subjects show significant M1VE-spine coherence\n', ...
    sum(sig_pos), nSubjects);

threshold      = 0.2;
prevalence_loc = mean(all_masks, 2);

% Subject line plot with ROI shading
subj_cmap = [27,158,119; 217,95,2; 117,112,179; 231,41,138; 102,166,30; ...
             230,171,2;  166,118,29] / 255;

hfig_lines = figure('Color','w'); hold on;

all_coh = cat(1, subjResults(:).coh_diff);
yl_pad  = [min(all_coh(:))*1.1, max(all_coh(:))*1.1];

fill([cord_pos(roi_idx(1))   cord_pos(roi_idx(end)) ...
      cord_pos(roi_idx(end)) cord_pos(roi_idx(1))], ...
     [yl_pad(1) yl_pad(1) yl_pad(2) yl_pad(2)], ...
     [0.85 0.85 0.85], 'EdgeColor','none', 'FaceAlpha',0.3, ...
     'HandleVisibility','off');

for s = 1:nSubjects
    cdiff = subjResults(s).coh_diff;
    thr   = subjResults(s).thr95;
    sig   = cdiff > thr;
    c     = subj_cmap(s,:);
    for i = 1:length(cord_pos)-1
        if sig(i) && sig(i+1)
            plot(cord_pos(i:i+1), cdiff(i:i+1), '-', 'Color',c, ...
                'LineWidth',1.5,'HandleVisibility','off');
        else
            plot(cord_pos(i:i+1), cdiff(i:i+1), '-', 'Color',[0.7 0.7 0.7], ...
                'LineWidth',1,'HandleVisibility','off');
        end
    end
    plot(cord_pos(sig), cdiff(sig), '.', 'Color',c, ...
        'MarkerSize',12,'HandleVisibility','off');
    if sig_pos(s)
        h(s) = plot(nan,nan,'-','Color',c,'LineWidth',1.5);
    else
        h(s) = plot(nan,nan,'-','Color',[0.7 0.7 0.7],'LineWidth',1.5);
    end
end
yline(0,':k','HandleVisibility','off');
ylim(yl_pad); xlim([min(cord_pos) max(cord_pos)]);
xlabel('Cranial-caudal position (mm)','FontSize',13);
ylabel('Coherence difference (stat-rest)','FontSize',13);
title('M1VE-spine coherence differences (true BEM, smoothed 20mm)');
legend(h, arrayfun(@(s) sprintf('Participant %d',s), 1:nSubjects,'UniformOutput',false), ...
    'Location','bestoutside');
set(gca,'FontSize',13); grid on;
savefig(hfig_lines, fullfile(fig_dir, ['stepB_group_lines' out_suffix '.fig']));
close(hfig_lines);

% Group prevalence line plot
hfig_prev = figure('Color','w','Position',[100 100 700 400]);
hold on;
yl_clust = [0, max(prevalence_loc)*1.1];
fill([cord_pos(roi_idx(1))   cord_pos(roi_idx(end)) ...
      cord_pos(roi_idx(end)) cord_pos(roi_idx(1))], ...
     [yl_clust(1) yl_clust(1) yl_clust(2) yl_clust(2)], ...
     [0.85 0.85 0.85], 'EdgeColor','none', 'FaceAlpha',0.3, ...
     'HandleVisibility','off');
plot(cord_pos, prevalence_loc, 'k-', 'LineWidth',2);
yline(threshold,'r--','LineWidth',1.5);
ylim(yl_clust); xlim([min(cord_pos) max(cord_pos)]);
xlabel('Cranial-caudal position (mm)','FontSize',13);
ylabel('Prevalence','FontSize',13);
title('Group prevalence — M1VE-spine coherence','Interpreter','none');
grid on; box on; set(gca,'FontSize',13);
savefig(hfig_prev, fullfile(fig_dir, ['stepB_group_prevalence' out_suffix '.fig']));
close(hfig_prev);

fprintf('\n=== STEP B DONE ===\n');

%% =========================================================================
%  LOCAL FUNCTIONS
%% =========================================================================
function [coh_diff, cohDiff_perm] = extract_coh_diff(source_perm, nsourcepoints, nPerm)
    cohDiff_perm = zeros(nsourcepoints, nPerm);
    for i = 1:nPerm
        cohDiff_perm(:,i) = source_perm.trialA(i).coh - source_perm.trialB(i).coh;
    end
    coh_diff = source_perm.avgA.coh - source_perm.avgB.coh;
end

function thr = compute_threshold(cohDiff_perm, mult_comp_corr, nsourcepoints)
    maxPerm = max(cohDiff_perm, [], 1);
    if mult_comp_corr
        thr = prctile(maxPerm, 95);
    else
        thr = prctile(cohDiff_perm, 95, 2);
        warning('Uncorrected per-source threshold used.');
    end
end

function pvals = compute_pvals(coh_diff, cohDiff_perm)
    nPerm = size(cohDiff_perm, 2);
    pvals = (sum(cohDiff_perm >= coh_diff, 2) + 1) / (nPerm + 1);
end

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