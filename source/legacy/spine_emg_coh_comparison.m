%% spine_emg_coherence_comparison.m

%%This script is a parameter comparison for the spine-EMG DICS coherence analysis for Participant 1 only, systematically varying three things:

%Leadfield model — BEM vs Biot-Savart (BSLaw), 2 options
%Regularisation (lambda) — 1%, 5%, 10%, 3 options
%Smoothing — on (20mm FWHM) vs off, 2 options
clear all; close all; clc;

%% =========================================================================
%  USER CONFIG
%% =========================================================================
fieldtrip_path = 'C:\Users\mspedden\Documents\fieldtrip';
spm_path       = 'C:\Users\mspedden\Documents\spm';
bsc_path       = 'C:\Users\mspedden\Documents\brainspineconnectivity\source';
data_root      = 'C:\spinecoh_data';
save_dir       = 'C:\Users\mspedden\Documents\brainspine_savetest\lf_comparison';

lf_configs(1).name   = 'BEM';
lf_configs(1).lf_path = 'C:\Leadfields meshes\leadfield_experimental_bem_experimental.mat';
lf_configs(1).lf_var  = 'leadfield_cord';

lf_configs(2).name   = 'BSLaw';
lf_configs(2).lf_path = 'C:\Leadfields meshes\leadfield_experimental_bslaw_experimental.mat';
lf_configs(2).lf_var  = 'leadfield_bs';

geomfile = 'C:\Leadfields meshes\geometries_experimental.mat';

smooth_vals = [0, 1];
lambda_vals = [1, 5, 10];
fwhm_mm     = 20;
radius_mm   = 3 * (fwhm_mm / 2.355);

sub            = 'OP00212';
fband          = [10 35];
numpermutation = 500;
roi_idx        = 25:30;
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

%% =========================================================================
%  LOAD GEOMETRY — shared source space and meshes
%% =========================================================================
fprintf('Loading geometry...\n');
geom_data  = load(geomfile);
sources_cent = geom_data.sources_cent;
mesh_torso   = geom_data.mesh_torso;

cord_pos      = sources_cent.pos(:,2);
nsourcepoints = size(sources_cent.pos, 1);
fprintf('  Source space: %d points, y range %.1f to %.1f mm\n', ...
    nsourcepoints, min(cord_pos), max(cord_pos));

%% =========================================================================
%  BUILD SMOOTHER ONCE (shared source space)
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
%  LOAD DATA ONCE
%% =========================================================================
fprintf('=== Loading and preprocessing data ===\n');
datafile = fullfile(data_root, ['sub-' sub], 'ses-001', 'meg', ...
    'pmergedoe1000mspddfflo45hi45hfcstatic_001_array1.mat');

D       = spm_eeg_load(datafile);
grad_mm = D.sensors('MEG');
ftdat   = spm2fieldtrip(D);

badchans = D.chanlabels(D.badchannels);
cfg = []; cfg.channel = setdiff(ftdat.label, badchans);
ftdat = ft_selectdata(cfg, ftdat);

cfg = []; cfg.rectify = 'yes'; cfg.channel = 'EXG1';
ftdatr = ft_preprocessing(cfg, ftdat);
for k = 1:length(ftdat.trial)
    ftdat.trial{k}(end,:) = ftdatr.trial{k};
end

trialinfo = ftdat.trialinfo;
statidx   = find(trialinfo == 1);
restidx   = find(trialinfo == 2);
nTrials   = min(numel(statidx), numel(restidx));

cfg = []; cfg.trials = statidx(1:nTrials);
statdat_tr = ft_selectdata(cfg, ftdat);
fprintf('  nTrials: %d\n', nTrials);

%% =========================================================================
%  VOLUME CONDUCTOR (shared)
%% =========================================================================
mesh_wm.unit = 'mm';
cfg = []; cfg.method = 'infinite'; cfg.siunits = 1;
cfg.grad = grad_mm; cfg.conductivity = 1;
dummyvol = ft_prepare_headmodel(cfg, mesh_torso);

%% =========================================================================
%  MAIN LOOP — leadfields x smoothing x regularisation
%% =========================================================================
n_lf     = numel(lf_configs);
n_smooth = numel(smooth_vals);
n_lambda = numel(lambda_vals);
n_total  = n_lf * n_smooth * n_lambda;

results  = struct();
cond_num = 0;

for li = 1:n_lf
    lfc = lf_configs(li);
    fprintf('\n=== Leadfield: %s ===\n', lfc.name);

    % Load leadfield
    lf_data = load(lfc.lf_path);
    lf_raw  = lf_data.(lfc.lf_var);

    % Label-based channel matching
    data_meg_labels        = ftdat.label(~strcmp(ftdat.label,'EXG1'));
    [common_labels,idx_lf] = intersect(lf_raw.label, data_meg_labels, 'stable');
    fprintf('  Data MEG: %d  |  LF: %d  |  Matched: %d\n', ...
        numel(data_meg_labels), numel(lf_raw.label), numel(common_labels));
    if numel(common_labels) < numel(data_meg_labels)
        fprintf('  WARNING: %d channels not in leadfield\n', ...
            numel(data_meg_labels) - numel(common_labels));
    end

    Lf       = lf_raw;
    Lf.label = common_labels;
    Lf.pos   = sources_cent.pos;   % shared source space in mm
    Lf.inside = ones(nsourcepoints, 1);
    for i = 1:numel(lf_raw.leadfield)
        if ~isempty(lf_raw.leadfield{i})
            Lf.leadfield{i} = lf_raw.leadfield{i}(idx_lf, :);
        end
    end

    % Frequency data
    cfg_sel = []; cfg_sel.channel = [Lf.label; {'EXG1'}];
    cfg_av  = []; cfg_av.avgoverfreq = 'yes';
    cfg_fr  = []; cfg_fr.output = 'powandcsd'; cfg_fr.method = 'mtmfft';
    cfg_fr.foilim = fband; cfg_fr.tapsmofrq = 1; cfg_fr.keeptrials = 'no';

    freqall = ft_freqanalysis(cfg_fr, ftdat);
    freqall = ft_selectdata(cfg_av,  freqall);
    freqall_lf = ft_selectdata(cfg_sel, freqall);

    cfg = []; cfg.trials = statidx(1:nTrials);
    statdat = ft_selectdata(cfg, ftdat);
    freqstat = ft_freqanalysis(cfg_fr, statdat);
    freqstat = ft_selectdata(cfg_av,  freqstat);
    freqstat_lf = ft_selectdata(cfg_sel, freqstat);

    cfg = []; cfg.trials = restidx(1:nTrials);
    restdat = ft_selectdata(cfg, ftdat);
    freqrest = ft_freqanalysis(cfg_fr, restdat);
    freqrest = ft_selectdata(cfg_av,  freqrest);
    freqrest_lf = ft_selectdata(cfg_sel, freqrest);

    statdat_shuf_base = statdat_tr;
    emg_chan_idx      = find(strcmp(statdat_tr.label, 'EXG1'));

    sourcemodel = [];
    sourcemodel.pos       = Lf.pos;
    sourcemodel.unit      = 'mm';
    sourcemodel.inside    = logical(Lf.inside);
    sourcemodel.leadfield = Lf.leadfield;
    sourcemodel.label     = Lf.label;

    % Inner loops
    for si = 1:n_smooth
        doSmooth = smooth_vals(si);

        for ri = 1:n_lambda
            lambda   = lambda_vals(ri);
            cond_num = cond_num + 1;

            cond_label = sprintf('%s_lam%d_smooth%d', lfc.name, lambda, doSmooth);
            fprintf('\n--- Condition %d/%d: %s ---\n', cond_num, n_total, cond_label);
            t_start = tic;

            % Common spatial filter
            cfg_dics = [];
            cfg_dics.sourcemodel     = sourcemodel;
            cfg_dics.headmodel       = dummyvol;
            cfg_dics.dics.keepfilter = 'yes';
            cfg_dics.dics.lambda     = sprintf('%d%%', lambda);
            cfg_dics.method          = 'dics';
            cfg_dics.refchan         = 'EXG1';
            coh_all = ft_sourceanalysis(cfg_dics, freqall_lf);

            cfg2 = cfg_dics;
            cfg2.dics.keepfilter = 'no';
            cfg2.dics.filter     = coh_all.avg.filter;

            source_stat = ft_sourceanalysis(cfg2, freqstat_lf);
            source_rest = ft_sourceanalysis(cfg2, freqrest_lf);

            coh_stat = source_stat.avg.coh;
            coh_rest = source_rest.avg.coh;
            coh_diff = coh_stat - coh_rest;

            if doSmooth
                coh_stat = Wsm * coh_stat;
                coh_rest = Wsm * coh_rest;
                coh_diff = Wsm * coh_diff;
            end

            % Permutation test
            fprintf('  Running %d permutations...\n', numpermutation);
            cohPerm = zeros(nsourcepoints, numpermutation);
            rng(1);

            for p = 1:numpermutation
                if mod(p,100)==0
                    fprintf('    Perm %d/%d\n', p, numpermutation);
                end
                emg_shuf     = randperm(nTrials);
                statdat_shuf = statdat_shuf_base;
                for tr = 1:nTrials
                    statdat_shuf.trial{tr}(emg_chan_idx,:) = ...
                        statdat_shuf_base.trial{emg_shuf(tr)}(emg_chan_idx,:);
                end

                cfg_fr2 = []; cfg_fr2.output = 'powandcsd'; cfg_fr2.method = 'mtmfft';
                cfg_fr2.foilim = fband; cfg_fr2.tapsmofrq = 1; cfg_fr2.keeptrials = 'no';
                freqshuf = ft_freqanalysis(cfg_fr2, statdat_shuf);
                cfg_av2  = []; cfg_av2.avgoverfreq = 'yes';
                freqshuf = ft_selectdata(cfg_av2, freqshuf);
                freqshuf = ft_selectdata(cfg_sel,  freqshuf);

                src_shuf = ft_sourceanalysis(cfg2, freqshuf);
                coh_p    = src_shuf.avg.coh;
                if doSmooth
                    coh_p = Wsm * coh_p;
                end
                cohPerm(:,p) = coh_p;
            end

            maxPerm = max(cohPerm, [], 1);
            thr95   = prctile(maxPerm, 95);
            mask    = coh_stat > thr95;

            fprintf('  Threshold (FWE p<0.05): %.6f\n', thr95);
            fprintf('  Significant sources:    %d / %d\n', sum(mask), nsourcepoints);
            [peak_coh, peak_idx] = max(coh_stat);
            fprintf('  Peak coherence: %.4f at y=%.1f mm (source %d)\n', ...
                peak_coh, cord_pos(peak_idx), peak_idx);
            fprintf('  Time: %.1f min\n', toc(t_start)/60);

            r = struct();
            r.cond_label  = cond_label;
            r.lf_name     = lfc.name;
            r.lambda      = lambda;
            r.doSmooth    = doSmooth;
            r.coh_stat    = coh_stat;
            r.coh_rest    = coh_rest;
            r.coh_diff    = coh_diff;
            r.cohPerm     = cohPerm;
            r.thr95       = thr95;
            r.mask        = mask;
            r.cord_pos    = cord_pos;
            results(cond_num) = r;

            save(fullfile(save_dir, ['result_' cond_label '.mat']), '-struct', 'r');
        end
    end
end

save(fullfile(save_dir, 'all_results.mat'), 'results', 'cord_pos', ...
    'nsourcepoints', 'roi_idx', 'fwhm_mm', 'lambda_vals', 'smooth_vals');
fprintf('\n\nAll conditions complete. Results saved.\n');

%% =========================================================================
%  FIGURES
%% =========================================================================
fprintf('Generating figures...\n');

lf_names   = {'BEM','BSLaw'};
cmap_lines = lines(3);

for si = 1:n_smooth
    doSmooth   = smooth_vals(si);
    smooth_str = {'No smoothing','Smoothed 20mm FWHM'};

    figure('Color','w','Position',[50 50 1400 700]);
    sgtitle(sprintf('Spine-EMG DICS coherence — %s', smooth_str{si}), ...
        'FontWeight','normal','FontSize',13);

    for li = 1:n_lf
        for ri = 1:n_lambda
            lambda = lambda_vals(ri);
            subplot_idx = (li-1)*n_lambda + ri;
            subplot(n_lf, n_lambda, subplot_idx);
            hold on;

            cond_label = sprintf('%s_lam%d_smooth%d', lf_names{li}, lambda, doSmooth);
            idx = find(strcmp({results.cond_label}, cond_label));
            if isempty(idx), continue; end
            r = results(idx);

            yl_pad = [min(r.coh_stat)*0.9, max(r.coh_stat)*1.15];
            fill([cord_pos(roi_idx(1))   cord_pos(roi_idx(end)) ...
                  cord_pos(roi_idx(end)) cord_pos(roi_idx(1))], ...
                 [yl_pad(1) yl_pad(1) yl_pad(2) yl_pad(2)], ...
                 [0.85 0.85 0.85], 'EdgeColor','none', 'DisplayName','ROI (C8-T1)');

            plot(cord_pos, r.coh_stat, '-', 'Color', cmap_lines(ri,:), ...
                'LineWidth', 2, 'DisplayName', 'Contraction');
            yline(r.thr95, '--', 'Color', cmap_lines(ri,:), 'LineWidth', 1.2, ...
                'DisplayName', sprintf('Thr (%.4f)', r.thr95));
            if any(r.mask)
                scatter(cord_pos(r.mask), r.coh_stat(r.mask), 40, ...
                    cmap_lines(ri,:), 'filled', 'DisplayName', 'Significant');
            end

            yline(0, 'k:', 'HandleVisibility','off');
            xlim([min(cord_pos) max(cord_pos)]);
            ylim(yl_pad);
            grid on; box on;

            title(sprintf('%s  \\lambda=%d%%', lf_names{li}, lambda), ...
                'FontWeight','normal','FontSize',10);
            if ri == 1, ylabel('Coherence','FontSize',9); end
            if li == n_lf
                xlabel('Position along cord (mm)','FontSize',9);
            else
                set(gca,'XTickLabel',[]);
            end
            if subplot_idx == 1
                legend('Location','northwest','FontSize',7);
            end
        end
    end

    fname = sprintf('comparison_smooth%d', doSmooth);
    savefig(gcf, fullfile(fig_dir, [fname '.fig']));
    saveas(gcf,  fullfile(fig_dir, [fname '.png']));
end

%% Summary bar chart
figure('Color','w','Position',[100 100 900 500]);
x_labels  = {};
peak_vals = zeros(1, n_total);
nsig_vals = zeros(1, n_total);

for k = 1:n_total
    r = results(k);
    [peak_vals(k), ~] = max(r.coh_stat);
    nsig_vals(k)      = sum(r.mask);
    x_labels{k}       = strrep(r.cond_label, '_', ' ');
end

subplot(1,2,1);
bar(peak_vals); set(gca,'XTick',1:n_total,'XTickLabel',x_labels,'FontSize',7);
xtickangle(35); ylabel('Peak coherence (contraction)');
title('Peak coherence','FontWeight','normal'); grid on; box on;

subplot(1,2,2);
bar(nsig_vals); set(gca,'XTick',1:n_total,'XTickLabel',x_labels,'FontSize',7);
xtickangle(35); ylabel('N significant sources');
title('Significant sources (FWE p<0.05)','FontWeight','normal'); grid on; box on;

sgtitle('Condition summary','FontWeight','normal','FontSize',12);
savefig(gcf, fullfile(fig_dir, 'summary_bar.fig'));
saveas(gcf,  fullfile(fig_dir, 'summary_bar.png'));
fprintf('Figures saved to %s\n', fig_dir);

%% =========================================================================
%  LOCAL FUNCTION
%% =========================================================================
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