%% spine_emg_coherence_bslaw.m
% Spine-EMG DICS coherence using BSLaw experimental leadfield
% Participant 1 only, with optional EMG trial shuffle permutation test.

clear all; close all; clc;

%% =========================================================================
%  USER CONFIG
%% =========================================================================
fieldtrip_path  = 'C:\Users\mspedden\Documents\fieldtrip';
spm_path        = 'C:\Users\mspedden\Documents\spm';
bsc_path        = 'C:\Users\mspedden\Documents\brainspineconnectivity\source';
data_root       = 'C:\spinecoh_data';
geomfile        = "C:\Leadfields meshes\geometries_experimental.mat";
lf_path         = 'C:\Leadfields meshes\leadfield_experimental_bslaw_experimental.mat';
save_dir        = 'C:\Users\mspedden\Documents\brainspine_savetest';

sub             = 'OP00212';
fband           = [10 35];
numpermutation  = 500;
rng(1);

% --- Flags ---
just_plot      = 0;   % Set 1 to skip analysis and plot from saved results
run_permtest   =1;   % Set 1 to run permutation test; 0 to just get peak coherence

%% =========================================================================
%  SETUP
%% =========================================================================
addpath(bsc_path); addpath(spm_path);
spm('defaults','EEG');
addpath(fieldtrip_path);
ft_defaults;

if ~exist(save_dir,'dir'), mkdir(save_dir); end

%% =========================================================================
%  LOAD GEOMETRY AND LEADFIELD
%% =========================================================================
fprintf('Loading geometry and leadfield...\n');
load(geomfile);
load(lf_path);

cord_pos      = sources_cent.pos(:,2);
nsourcepoints = size(sources_cent.pos, 1);

if just_plot
    load(fullfile(save_dir, 'spine_emg_coherence_bslaw_p1.mat'));
    load(geomfile);
    nsourcepoints = size(sources_cent.pos, 1);
    cord_pos      = sources_cent.pos(:,2);
else

%% =========================================================================
%  LOAD AND PREPROCESS DATA
%% =========================================================================
fprintf('Loading data for %s...\n', sub);
datafile = fullfile(data_root, ['sub-' sub], 'ses-001', 'meg', ...
    'pmergedoe1000mspddfflo45hi45hfcstatic_001_array1.mat');

D       = spm_eeg_load(datafile);
grad_mm = D.sensors('MEG');
ftdat   = spm2fieldtrip(D);

% Remove bad channels
badchans = D.chanlabels(D.badchannels);
cfg = []; cfg.channel = setdiff(ftdat.label, badchans);
ftdat = ft_selectdata(cfg, ftdat);

% Rectify EMG
cfg = []; cfg.rectify = 'yes'; cfg.channel = 'EXG1';
ftdatr = ft_preprocessing(cfg, ftdat);
for k = 1:length(ftdat.trial)
    ftdat.trial{k}(end,:) = ftdatr.trial{k};
end

%% =========================================================================
%  CHANNEL MATCHING — label-based (BSLaw leadfield labels match data)
%% =========================================================================
fprintf('Matching channels by label...\n');

% Build a leadfield struct with only channels present in the data
% (minus bad channels), matched by label
data_meg_labels = ftdat.label(~strcmp(ftdat.label, 'EXG1'));  % MEG only
lf_labels       = leadfield_bs.label;

% Find the subset of leadfield rows that correspond to data channels
[common_labels, idx_lf, idx_dat] = intersect(lf_labels, data_meg_labels, 'stable');

fprintf('  Data MEG channels:      %d\n', numel(data_meg_labels));
fprintf('  Leadfield channels:     %d\n', numel(lf_labels));
fprintf('  Matched (intersection): %d\n', numel(common_labels));

if numel(common_labels) < numel(data_meg_labels)
    missing = setdiff(data_meg_labels, lf_labels);
    fprintf('  WARNING: %d data channels not in leadfield:\n', numel(missing));
    fprintf('    %s\n', missing{:});
end

% Build trimmed leadfield using matched rows only
Lf       = leadfield_bs;
Lf.label = common_labels;
for i = 1:numel(leadfield_bs.leadfield)
    if ~isempty(leadfield_bs.leadfield{i})
        Lf.leadfield{i} = leadfield_bs.leadfield{i}(idx_lf, :);
    end
end
fprintf('  Leadfield built: %d channels x sources\n', numel(Lf.label));

%% =========================================================================
%  VOLUME CONDUCTOR
%% =========================================================================
mesh_wm.unit = 'mm';
cfg = []; cfg.method = 'infinite'; cfg.siunits = 1;
cfg.grad = grad_mm; cfg.conductivity = 1;
dummyvol = ft_prepare_headmodel(cfg, mesh_torso);

%% =========================================================================
%  FREQUENCY ANALYSIS
%% =========================================================================
fprintf('Computing frequency data...\n');
trialinfo = ftdat.trialinfo;
statidx   = find(trialinfo == 1);
restidx   = find(trialinfo == 2);
nTrials   = min(numel(statidx), numel(restidx));

% Channel selection config — MEG channels + EMG, reused throughout
cfg_sel = []; cfg_sel.channel = [Lf.label; {'EXG1'}];

cfg_av = []; cfg_av.avgoverfreq = 'yes';  % reused below

% Combined — for common spatial filter
cfg = []; cfg.output = 'powandcsd'; cfg.method = 'mtmfft';
cfg.foilim = fband; cfg.tapsmofrq = 1; cfg.keeptrials = 'no';
freqall     = ft_freqanalysis(cfg, ftdat);
freqall     = ft_selectdata(cfg_av,  freqall);
freqall_lf  = ft_selectdata(cfg_sel, freqall);

% Contraction
cfg = []; cfg.trials = statidx(1:nTrials);
statdat = ft_selectdata(cfg, ftdat);
cfg = []; cfg.output = 'powandcsd'; cfg.method = 'mtmfft';
cfg.foilim = fband; cfg.tapsmofrq = 1; cfg.keeptrials = 'no';
freqstat    = ft_freqanalysis(cfg, statdat);
freqstat    = ft_selectdata(cfg_av,  freqstat);
freqstat_lf = ft_selectdata(cfg_sel, freqstat);

% Rest
cfg = []; cfg.trials = restidx(1:nTrials);
restdat = ft_selectdata(cfg, ftdat);
cfg = []; cfg.output = 'powandcsd'; cfg.method = 'mtmfft';
cfg.foilim = fband; cfg.tapsmofrq = 1; cfg.keeptrials = 'no';
freqrest    = ft_freqanalysis(cfg, restdat);
freqrest    = ft_selectdata(cfg_av,  freqrest);
freqrest_lf = ft_selectdata(cfg_sel, freqrest);

Lf.inside=ones(length(Lf.pos), 1);

%% =========================================================================
%  DICS — COMMON SPATIAL FILTER FROM COMBINED DATA
%% =========================================================================
fprintf('Computing common spatial filter...\n');
cfg = [];
cfg.sourcemodel.pos       = Lf.pos;
cfg.sourcemodel.unit      = 'mm';
cfg.sourcemodel.inside    = logical(Lf.inside);
cfg.sourcemodel.leadfield = Lf.leadfield;
cfg.sourcemodel.label     = Lf.label;
cfg.headmodel             = dummyvol;
cfg.dics.keepfilter       = 'yes';
cfg.dics.lambda           = 10;
cfg.method                = 'dics';
cfg.refchan               = 'EXG1';
coh_all = ft_sourceanalysis(cfg, freqall_lf);

% Config for applying fixed filter
cfg2 = [];
cfg2.sourcemodel.pos       = Lf.pos;
cfg2.sourcemodel.unit      = 'mm';
cfg2.sourcemodel.inside    = logical(Lf.inside);
cfg2.sourcemodel.leadfield = Lf.leadfield;
cfg2.sourcemodel.label     = Lf.label;
cfg2.headmodel             = dummyvol;
cfg2.dics.filter           = coh_all.avg.filter;
cfg2.dics.lambda           = 10;
cfg2.method                = 'dics';
cfg2.refchan               = 'EXG1';

% Apply to contraction and rest
source_stat = ft_sourceanalysis(cfg2, freqstat_lf);
source_rest = ft_sourceanalysis(cfg2, freqrest_lf);

coh_stat = source_stat.avg.coh;
coh_rest = source_rest.avg.coh;
coh_diff = coh_stat - coh_rest;

[peak_coh, peak_idx] = max(coh_stat);
fprintf('\n--- Peak coherence (no permutation test) ---\n');
fprintf('  Peak coherence (contraction): %.4f\n',  peak_coh);
fprintf('  Peak position:                y=%.1f mm (source %d)\n', ...
    cord_pos(peak_idx), peak_idx);
fprintf('  Max coherence (rest):         %.4f\n',  max(coh_rest));
fprintf('  Max coherence (diff):         %.4f\n\n', max(coh_diff));

%% =========================================================================
%  OPTIONAL PERMUTATION TEST
%% =========================================================================
if run_permtest
    fprintf('Running %d permutations (EMG trial shuffle)...\n', numpermutation);

    cfg = []; cfg.trials = statidx(1:nTrials);
    statdat_tr = ft_selectdata(cfg, ftdat);

    emg_chan_idx = find(strcmp(statdat_tr.label, 'EXG1'));
    cohPerm      = zeros(nsourcepoints, numpermutation);

    for p = 1:numpermutation
        if mod(p,50)==0, fprintf('  Permutation %d/%d\n', p, numpermutation); end

        emg_shuf     = randperm(nTrials);
        statdat_shuf = statdat_tr;
        for tr = 1:nTrials
            statdat_shuf.trial{tr}(emg_chan_idx,:) = ...
                statdat_tr.trial{emg_shuf(tr)}(emg_chan_idx,:);
        end

        cfg_fr = []; cfg_fr.output = 'powandcsd'; cfg_fr.method = 'mtmfft';
        cfg_fr.foilim = fband; cfg_fr.tapsmofrq = 1; cfg_fr.keeptrials = 'no';
        freqshuf = ft_freqanalysis(cfg_fr, statdat_shuf);
        cfg_av2  = []; cfg_av2.avgoverfreq = 'yes';
        freqshuf = ft_selectdata(cfg_av2, freqshuf);
        freqshuf = ft_selectdata(cfg_sel, freqshuf);

        src_shuf     = ft_sourceanalysis(cfg2, freqshuf);
        cohPerm(:,p) = src_shuf.avg.coh;
    end

    maxPerm = max(cohPerm, [], 1);
    thr95   = prctile(maxPerm, 95);
    mask    = coh_stat > thr95;

    fprintf('\n--- Permutation test results ---\n');
    fprintf('  Threshold (FWE p<0.05): %.6f\n', thr95);
    fprintf('  Significant sources:    %d / %d\n', sum(mask), nsourcepoints);
    if any(mask)
        [~, peak_idx] = max(coh_stat .* double(mask));
        fprintf('  Peak position:          y=%.1f mm (source %d)\n', ...
            cord_pos(peak_idx), peak_idx);
    end
else
    % Dummy values so save/plot code works whether or not permtest ran
    cohPerm = [];
    thr95   = NaN;
    mask    = false(nsourcepoints, 1);
    fprintf('Permutation test skipped (run_permtest=0).\n');
end

%% Save
save(fullfile(save_dir, 'spine_emg_coherence_bslaw_p1.mat'), ...
    'coh_stat', 'coh_rest', 'coh_diff', 'cohPerm', ...
    'thr95', 'mask', 'cord_pos');
fprintf('Saved.\n');
end  % end just_plot

%% =========================================================================
%  FIGURES
%% =========================================================================

% --- Line plot ---
figure('Color','w','Position',[100 100 750 500]);
hold on;

% ROI shading — C8-T1 cord segments (vertebral levels C6-T1), source indices 25-30
roi_idx = 25:30;
yl_pad  = [min(coh_stat)*0.9, max(coh_stat)*1.15];
fill([cord_pos(roi_idx(1))   cord_pos(roi_idx(end)) ...
      cord_pos(roi_idx(end)) cord_pos(roi_idx(1))], ...
     [yl_pad(1) yl_pad(1) yl_pad(2) yl_pad(2)], ...
     [0.85 0.85 0.85], 'EdgeColor', 'none', 'DisplayName', 'ROI (C8-T1)');

plot(cord_pos, coh_stat, 'b-', 'LineWidth', 2, 'DisplayName', 'Contraction');
if run_permtest && ~isnan(thr95)
    yline(thr95, 'b--', 'LineWidth', 1.5, ...
        'DisplayName', sprintf('Threshold (%.4f)', thr95));
    scatter(cord_pos(mask), coh_stat(mask), 60, 'b', 'filled', ...
        'DisplayName', 'Significant');
end
yline(0, 'k:', 'HandleVisibility','off');
ylabel('Coherence');
xlabel('Position along cord (mm, inferior \rightarrow superior)');
title('Participant 1 — Spine-EMG coherence (BSLaw)', 'FontWeight','normal');
legend('Location','northwest','FontSize',10);
grid on; box on;
xlim([min(cord_pos) max(cord_pos)]);
ylim(yl_pad);

% --- Null maxima histogram (only if permutation test was run) ---
if run_permtest && ~isempty(cohPerm)
    [~, maxIdx_perm] = max(cohPerm, [], 1);
    [~, obsMaxIdx]   = max(coh_stat);

    figure('Color','w','Position',[100 100 600 450]);
    hold on;
    histogram(cord_pos(maxIdx_perm), 44, ...
        'FaceColor',[0.75 0.75 0.75],'EdgeColor','k','LineWidth',0.8);
    xline(cord_pos(obsMaxIdx), '-', 'Color',[0.2 0 0], 'LineWidth', 2);
    xlabel('Cranio-caudal position (mm)','FontSize',14);
    ylabel('Count','FontSize',14);
    legend({'Null maxima','Observed maximum'},'Location','best','FontSize',12,'Box','off');
    set(gca,'FontSize',14,'LineWidth',1.2,'TickDir','out'); box off;
    title('OP00212 — null maxima (EMG trial shuffle)','Interpreter','none');
end

% --- Coherence on bone mesh ---
coh_norm = (coh_stat - min(coh_stat)) / (max(coh_stat) - min(coh_stat));
cmap_jet = jet(256);
colors   = interp1(linspace(0,1,256), cmap_jet, coh_norm);

figure('Color','w','Position',[100 100 700 700]);
ft_plot_mesh(mesh_bone, 'facecolor', [0.9 0.85 0.7], ...
    'facealpha', 0.3, 'edgecolor', 'none');
hold on;
for s = 1:nsourcepoints
    plot3(sources_cent.pos(s,1), sources_cent.pos(s,2), sources_cent.pos(s,3), ...
        'o', 'MarkerFaceColor', colors(s,:), 'MarkerEdgeColor', 'k', ...
        'MarkerSize', 12, 'LineWidth', 0.5);
end
for s = 1:nsourcepoints
    text(sources_cent.pos(s,1), sources_cent.pos(s,2), ...
        sources_cent.pos(s,3)+2, num2str(s), ...
        'FontSize', 9, 'Color', 'k', 'HorizontalAlignment', 'center');
end
colormap(jet); clim([min(coh_stat) max(coh_stat)]);
cb = colorbar; cb.Label.String = 'Coherence (contraction)';
view(-250,-1); camlight; lighting gouraud;
title('Participant 1 — contraction coherence on vertebrae (BSLaw)', ...
    'FontWeight','normal','FontSize',12);