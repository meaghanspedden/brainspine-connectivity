%% brain_lambda_10pct_all.m
% Brain-EMG DICS coherence with lambda='10%' for all 7 brain subjects.
% Corrected regularisation matching spine analysis.

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

subs_brain = {'OP00212','OP00213','OP00215','OP00219', ...
    'OP00225','OP00221','OP00224'};

% --- Flags ---
run_subjects   = 1;   % set 1 to run subject loop
run_group_only = 0;   % set 1 to just load results and run group

% --- Analysis settings --

fband          = [10 35];
numpermutation = 500;
lambda         = '10%';
fwhm_mm        = 8;
radius_mm      = 3 * (fwhm_mm / 2.355);
out_suffix     = '_brain_pct10';
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

invpthr = -log10(0.05);

%% =========================================================================
%  LOAD GEOMETRY ONCE
%% =========================================================================
fprintf('Loading geometry...\n');
load(geomfile);   % sources_brain, mesh_brain, mesh_torso etc.
mesh_brain.unit = 'mm';
load(T_mat_path); T_inv = inv(T);

%% =========================================================================
%  SUBJECT LOOP
%% =========================================================================
subjResults = struct();
if run_subjects
    for ss = 1:7
        sub = subs_brain{ss};
        fprintf('\n========================================\n');
        fprintf('  Subject %s (%d/%d)\n', sub, ss, length(subs_brain));
        fprintf('========================================\n');
        t_sub = tic;

        %% Load data
        run = '001'; if strcmp(sub,'OP00224'), run = '002'; end
        datafile = fullfile(data_root, ['sub-' sub], 'ses-001', 'meg', ...
            sprintf('pmergedoe1000mspddfflo45hi45hfcstatic_%s_array1.mat', run));

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

        %% Brain sensors and forward model
%         brainidx  = find(grad_mm.chanpos(:,2) > 200);
%         braingrad = subset_grad(grad_mm, brainidx);

        cfg = []; cfg.method = 'singleshell';
        vol = ft_prepare_headmodel(cfg, mesh_brain);
        cfg = []; cfg.sourcemodel = sources_brain; cfg.headmodel = vol;
        cfg.grad = grad_mm; cfg.reducerank = 'no';
        LF = ft_prepare_leadfield(cfg);

        %% Frequency data
        [statdat, restdat, freqdat] = make_freq_data(ftdat, ftdat.trialinfo, fband);

        %% Common spatial filter
        fprintf('  Computing common spatial filter...\n');
        cfg_dics = [];
        cfg_dics.grid                  = sources_brain;
        cfg_dics.headmodel             = vol;
        cfg_dics.sourcemodel.leadfield = LF.leadfield;
        cfg_dics.dics.keepfilter       = 'yes';
        cfg_dics.dics.lambda           = lambda;
        cfg_dics.method                = 'dics';
        cfg_dics.refchan               = 'EXG1';
        coh_source = ft_sourceanalysis(cfg_dics, freqdat);

        %% Permutation test
        fprintf('  Running permutation test (%d permutations)...\n', numpermutation);
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

        %% Build smoother on first subject
        if ss == 1
            fprintf('Building brain smoother (FWHM=%d mm)...\n', fwhm_mm);
            inside_idx    = (1:size(source_perm.pos,1))';
            inside_pos    = source_perm.pos;
            Wsm           = make_gaussian_smoother(inside_pos, fwhm_mm, radius_mm);

        end
        nsourcepoints = size(source_perm.pos,1);
        fprintf('  Sources: %d\n', nsourcepoints);
        nPerm = numel(source_perm.trialA);
        [coh_diff, cohDiff_perm] = extract_coh_diff_brain(source_perm, nsourcepoints, nPerm);

        %% Smooth
        coh_diff_inside            = coh_diff(inside_idx);
        cohDiff_perm_inside        = cohDiff_perm(inside_idx,:);
        coh_diff_inside            = Wsm * coh_diff_inside;
        cohDiff_perm_inside        = Wsm * cohDiff_perm_inside;
        coh_diff(inside_idx)       = coh_diff_inside;
        cohDiff_perm(inside_idx,:) = cohDiff_perm_inside;

        %% Threshold
        maxPerm = max(cohDiff_perm_inside, [], 1);
        thr95   = prctile(maxPerm, 95);
        mask    = coh_diff > thr95;

        pvals       = (sum(cohDiff_perm >= coh_diff, 2) + 1) / (nPerm + 1);
        invp        = -log10(pvals);
        invp_masked = invp; invp_masked(~mask) = 0;

        fprintf('  Threshold: %.6f  |  Sig sources: %d / %d\n', ...
            thr95, sum(mask), nsourcepoints);
        fprintf('  Subject time: %.1f min\n', toc(t_sub)/60);

        %% MNI peak
        maxval   = max(coh_diff);
        maxidx   = find(coh_diff == maxval);
        maxpos   = source_perm.pos(maxidx,:);
        maxpos_h = [maxpos, ones(size(maxpos,1),1)]';
        x_mni    = (T_inv * maxpos_h)'; x_mni = x_mni(:,1:3);
        fprintf('  Peak MNI: [%.1f %.1f %.1f]\n', x_mni(1,:));

        %% Interpolate onto brain mesh
        source_p         = coh_source;
        source_p.avg.coh = invp_masked;
        cfg_interp = []; cfg_interp.parameter = 'coh';
        cfg_interp.interpmethod = 'sphere_avg';
        cfg_interp.sphereradius = 10;
        brain_int = ft_sourceinterpolate(cfg_interp, source_p, mesh_brain);

        %% Figure — individual subject
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
        title(sprintf('%s — brain-EMG coherence (lambda=10%%)', sub),'Interpreter','none');
        savefig(hfig, fullfile(fig_dir, ...
            sprintf('step1_sub%s_brainEMG_coherence%s.fig', sub, out_suffix)));
        close(hfig);

        %% Store
        subjResults(ss).subjID    = sub;
        subjResults(ss).coh_diff  = coh_diff;
        subjResults(ss).thr95     = thr95;
        subjResults(ss).sig_mask  = double(mask);
        subjResults(ss).pos       = source_perm.pos;
        subjResults(ss).inside    = true(nsourcepoints,1);
        subjResults(ss).invp      = invp;
        subjResults(ss).x_mni     = x_mni;
        subjResults(ss).brain_int = brain_int;

        save(fullfile(save_dir, sprintf('brainResult_sub%s%s.mat', sub, out_suffix)), ...
            'coh_diff','cohDiff_perm','thr95','mask','invp','x_mni');
    end
    %% Save group results
    save(fullfile(save_dir, ['groupRes_brain_DICS' out_suffix '.mat']), 'subjResults');
    fprintf('\nAll subjects complete. \n');
end


if run_group_only
    fprintf('Loading saved results and running group analysis...\n');
    loaded = load(fullfile(save_dir, ['groupRes_brain_DICS' out_suffix '.mat']));
    subjResults = loaded.subjResults;
    % Add pos and inside — same grid for all subjects
    for ss = 1:numel(subjResults)
        subjResults(ss).pos    = sources_brain;
        subjResults(ss).inside = true(size(subjResults(ss).coh_diff));
    end

end

%% =========================================================================
%  GROUP ANALYSIS
%% =========================================================================
run_group_brain(subjResults, geomfile, T_mat_path, cmap_hot, fwhm_mm, ...
    fig_dir, out_suffix);

fprintf('\n=== DONE ===\n');

%% =========================================================================
%  LOCAL FUNCTIONS
%% =========================================================================

function run_group_brain(subjResults, geomfile, T_mat_path, cmap, fwhm, ...
    fig_dir, out_suffix)

load(geomfile); mesh_brain.unit = 'mm';
load(T_mat_path); T_inv = inv(T);

nSubs     = length(subjResults);
all_masks = cat(2, subjResults(:).sig_mask);

sig_pos = false(nSubs,1);
for ss = 1:nSubs
    if any(subjResults(ss).coh_diff > subjResults(ss).thr95)
        sig_pos(ss) = true;
    end
end
fprintf('  %d/%d subjects show significant brain-EMG coherence (smoothed %d mm, lambda=10%%)\n', ...
    sum(sig_pos), nSubs, fwhm);

group_prevalence = mean(all_masks, 2);
threshold        = 0.3;
gp_masked        = group_prevalence;
gp_masked(gp_masked < threshold) = 0;

group_ft         = subjResults(1).brain_int;
group_ft.pow     = gp_masked;
group_ft.pow     = gp_masked;

% Group MNI peak
maxval   = max(group_ft.pow);
maxidx   = find(group_ft.pow == maxval);
maxpos   = group_ft.pos(maxidx,:);
maxpos_h = [maxpos, ones(size(maxpos,1),1)]';
x_mni    = (T_inv * maxpos_h)'; x_mni = x_mni(:,1:3);
fprintf('  Group peak MNI: [%.1f %.1f %.1f]\n', x_mni(1,:));
if length(maxidx) > 1, disp('  Multiple group max locations'); end

% Group max location on mesh
hfig_max = figure;
ft_plot_mesh(mesh_brain); hold on;
plot3(maxpos(:,1), maxpos(:,2), maxpos(:,3), 'r*', 'MarkerSize', 12);
title(sprintf('Group max — brain-EMG coherence (lambda=10%%)'), 'Interpreter','none');
savefig(hfig_max, fullfile(fig_dir, ['step1_group_brain_max' out_suffix '.fig']));
close(hfig_max);

% Interpolate
cfg = []; cfg.parameter = 'pow';
cfg.interpmethod = 'sphere_avg'; cfg.sphereradius = 10;
group_int = ft_sourceinterpolate(cfg, group_ft, mesh_brain);

% Group prevalence figure
ncol     = 256;
hfig_prev = figure;
cfg2 = []; cfg2.method = 'surface'; cfg2.funparameter = 'pow';
cfg2.funcolorlim = [threshold max(group_int.pow)];
cfg2.funcolormap = cmap;
cfg2.projmethod  = 'nearest';
cfg2.surffile    = mesh_brain;
ft_sourceplot(cfg2, group_int);
view(176,-10); camlight; ax = gca; ax.FontSize = 14;
hpatch = findobj(gcf,'Type','patch'); set(hpatch,'FaceAlpha',0.9);
title(sprintf('Group prevalence — brain-EMG coherence (smoothed %d mm, lambda=10%%)', fwhm), ...
    'Interpreter','none');
savefig(hfig_prev, fullfile(fig_dir, ['step1_group_brain_prevalence' out_suffix '.fig']));


end

% -------------------------------------------------------------------------
function [coh_diff, cohDiff_perm] = extract_coh_diff_brain(source_perm, nsourcepoints, nPerm)
cohDiff_perm = zeros(nsourcepoints, nPerm);
for i = 1:nPerm
    cohDiff_perm(:,i) = source_perm.trialA(i).coh - source_perm.trialB(i).coh;
end
coh_diff = source_perm.avgA.coh - source_perm.avgB.coh;
end

% -------------------------------------------------------------------------
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

% -------------------------------------------------------------------------
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