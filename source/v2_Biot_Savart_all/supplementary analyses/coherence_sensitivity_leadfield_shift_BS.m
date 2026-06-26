%% coherence_shift_sensitivity_sub1.m
% Runs spine-EMG DICS coherence for participant 1 with original and one
% randomly selected ~15mm shifted BS leadfield, then compares peak
% coherence location along the cord.

clear all; close all; clc;

%% =========================================================================
%  USER CONFIG
%% =========================================================================
fieldtrip_path = 'C:\Users\mspedden\Documents\fieldtrip';
spm_path       = 'C:\Users\mspedden\Documents\spm';
bsc_path       = 'C:\Users\mspedden\Documents\brainspineconnectivity\source';
data_root      = 'C:\spinecoh_data';
save_dir       = 'C:\Users\mspedden\Documents\brainspine_savetest';

geomfile       = 'C:\Leadfields meshes\geometries_experimental_withbrain.mat';
lf_path_orig   = 'C:\Leadfields meshes\leadfield_experimental_bslaw_experimental.mat';
bs_lf_path     = 'C:\Users\mspedden\Documents\bslaw_sensitivity_analysis\bslaw_sensitivity_analysis\bs_law_fields';

sub            = 'OP00212';
fband          = [10 35];
numpermutation = 500;
lambda         = '10%';
fwhm_mm        = 20;
radius_mm      = 3 * (fwhm_mm / 2.355);

%% =========================================================================
%  RANDOMLY SELECT ONE ~15mm SHIFT CONDITION
%% =========================================================================
shift_candidates = { ...
    'sensor_bundle3_shift1',  15.67; ...
    'sensor_bundle3_shift3',  14.22; ...
    'sensor_bundle3_shift6',  17.29; ...
    'sensor_bundle3_shift7',  17.93  };

rng('shuffle');
chosen_idx   = randi(size(shift_candidates, 1));
chosen_label = shift_candidates{chosen_idx, 1};
chosen_mag   = shift_candidates{chosen_idx, 2};
fprintf('Selected shift condition: %s (|shift|=%.2f mm)\n', chosen_label, chosen_mag);

lf_path_shift = fullfile(bs_lf_path, ...
    sprintf('leadfield_%s_bslaw_experimental.mat', chosen_label));
assert(exist(lf_path_shift,'file')==2, 'Shifted leadfield not found: %s', lf_path_shift);

%% =========================================================================
%  SETUP
%% =========================================================================
addpath(bsc_path); addpath(spm_path);
spm('defaults','EEG');
addpath(fieldtrip_path);
ft_defaults;

fig_dir = fullfile(save_dir, 'figures');
if ~exist(fig_dir,'dir'), mkdir(fig_dir); end

%% =========================================================================
%  LOAD GEOMETRY
%% =========================================================================
fprintf('Loading geometry...\n');
geom_exp     = load(geomfile);
sources_cent = geom_exp.sources_cent;
mesh_torso   = geom_exp.mesh_torso;

nsourcepoints = size(sources_cent.pos, 1);

%% =========================================================================
%  BUILD SMOOTHER
%% =========================================================================
fprintf('Building Gaussian smoother (FWHM=%d mm)...\n', fwhm_mm);
Wsm = make_gaussian_smoother(sources_cent.pos, fwhm_mm, radius_mm);

%% =========================================================================
%  LOAD SUBJECT DATA
%% =========================================================================
fprintf('Loading subject %s...\n', sub);
datafile = fullfile(data_root, ['sub-' sub], 'ses-001', 'meg', ...
    sprintf('pmergedoe1000mspddfflo45hi45hfcstatic_001_array1.mat'));

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

% Trial indices
trialinfo = ftdat.trialinfo;
statidx   = find(trialinfo == 1);
restidx   = find(trialinfo == 2);
nTrials   = min(numel(statidx), numel(restidx));
fprintf('  nTrials (stat/rest): %d\n', nTrials);

% Frequency decomposition — done once, shared across both runs
cfg_fr = []; cfg_fr.output = 'powandcsd'; cfg_fr.method = 'mtmfft';
cfg_fr.foilim = fband; cfg_fr.tapsmofrq = 1; cfg_fr.keeptrials = 'yes';
cfg_av = []; cfg_av.avgoverfreq = 'yes';
cfg_sel_emg = []; cfg_sel_emg.channel = [ftdat.label(~strcmp(ftdat.label,'EXG1')); {'EXG1'}];

freqdat_tr = ft_freqanalysis(cfg_fr, ftdat);
freqdat_tr = ft_selectdata(cfg_av, freqdat_tr);

cfg = []; cfg.trials = statidx(1:nTrials);
statdat_full = ft_selectdata(cfg, freqdat_tr);
cfg = []; cfg.trials = restidx(1:nTrials);
restdat_full = ft_selectdata(cfg, freqdat_tr);

cfg_fr2 = []; cfg_fr2.output = 'powandcsd'; cfg_fr2.method = 'mtmfft';
cfg_fr2.foilim = fband; cfg_fr2.tapsmofrq = 1; cfg_fr2.keeptrials = 'no';
cfg_av2 = []; cfg_av2.avgoverfreq = 'yes';
freqdat_comb = ft_freqanalysis(cfg_fr2, ftdat);
freqdat_comb = ft_selectdata(cfg_av2, freqdat_comb);

% Volume conductor (infinite medium — same for both runs)
cfg_vol = []; cfg_vol.method = 'infinite'; cfg_vol.siunits = 1;
cfg_vol.grad = grad_mm; cfg_vol.conductivity = 1;
dummyvol = ft_prepare_headmodel(cfg_vol, mesh_torso);

%% =========================================================================
%  RUN COHERENCE — ORIGINAL AND SHIFTED
%% =========================================================================
lf_paths    = {lf_path_orig,  lf_path_shift};
run_labels  = {'Original',    sprintf('Shifted (%.1f mm)', chosen_mag)};
results_coh = struct();

for run = 1:2
    fprintf('\n── %s leadfield ──\n', run_labels{run});

    lf_data = load(lf_paths{run});
    lf_raw  = lf_data.leadfield_bs;

    data_meg_labels        = ftdat.label(~strcmp(ftdat.label,'EXG1'));
    [common_labels,idx_lf] = intersect(lf_raw.label, data_meg_labels, 'stable');
    fprintf('  Matched channels: %d\n', numel(common_labels));

    Lf          = lf_raw;
    Lf.label    = common_labels;
    Lf.pos      = sources_cent.pos;
    Lf.inside   = ones(nsourcepoints, 1);
    for i = 1:numel(lf_raw.leadfield)
        if ~isempty(lf_raw.leadfield{i})
            Lf.leadfield{i} = lf_raw.leadfield{i}(idx_lf, :);
        end
    end

    cfg_sel = []; cfg_sel.channel = [Lf.label; {'EXG1'}];
    statdat  = ft_selectdata(cfg_sel, statdat_full);
    restdat  = ft_selectdata(cfg_sel, restdat_full);
    freqdat  = ft_selectdata(cfg_sel, freqdat_comb);

    sourcemodel           = [];
    sourcemodel.pos       = Lf.pos;
    sourcemodel.unit      = 'mm';
    sourcemodel.inside    = logical(Lf.inside);
    sourcemodel.leadfield = Lf.leadfield;
    sourcemodel.label     = Lf.label;

    % Common spatial filter
    cfg_dics = [];
    cfg_dics.sourcemodel     = sourcemodel;
    cfg_dics.headmodel       = dummyvol;
    cfg_dics.dics.keepfilter = 'yes';
    cfg_dics.dics.lambda     = lambda;
    cfg_dics.method          = 'dics';
    cfg_dics.refchan         = 'EXG1';
    coh_source = ft_sourceanalysis(cfg_dics, freqdat);

    % Permutation test
    fprintf('  Running %d permutations...\n', numpermutation);
    cfg_perm = [];
    cfg_perm.sourcemodel        = sourcemodel;
    cfg_perm.headmodel          = dummyvol;
    cfg_perm.dics.filter        = coh_source.avg.filter;
    cfg_perm.dics.lambda        = lambda;
    cfg_perm.method             = 'dics';
    cfg_perm.refchan            = 'EXG1';
    cfg_perm.permutation        = 'yes';
    cfg_perm.numpermutation     = numpermutation;
    source_perm = ft_sourceanalysis(cfg_perm, statdat, restdat);

    nPerm = numel(source_perm.trialA);
    [coh_diff, cohDiff_perm] = extract_coh_diff(source_perm, nsourcepoints, nPerm);

    cohDiff_perm = Wsm * cohDiff_perm;
    coh_diff     = Wsm * coh_diff;

    thr95 = compute_threshold(cohDiff_perm, nsourcepoints);
    mask  = coh_diff > thr95;

    invp_smooth  = smooth_invp(coh_diff, cohDiff_perm, nsourcepoints, nPerm);
    [peak_val, peak_idx] = max(invp_smooth);
    peak_pos = sources_cent.pos(peak_idx, 2);

    fprintf('  Threshold: %.6f\n', thr95);
    fprintf('  Significant sources: %d / %d\n', sum(mask), nsourcepoints);
    fprintf('  Peak -log10(p): %.3f at y=%.1f mm\n', peak_val, peak_pos);

    results_coh(run).label       = run_labels{run};
    results_coh(run).coh_diff    = coh_diff;
    results_coh(run).invp_smooth = invp_smooth;
    results_coh(run).thr95       = thr95;
    results_coh(run).mask        = mask;
    results_coh(run).peak_pos    = peak_pos;
    results_coh(run).peak_val    = peak_val;
    results_coh(run).peak_idx    = peak_idx;
end

%% =========================================================================
%  FIGURE
%% =========================================================================
x      = sources_cent.pos(:, 2);
cols   = [0.02 0.15 0.50;   % dark blue  — original
          0.80 0.20 0.10];  % red-orange — shifted

figure('Color','w','Position',[100 100 700 420]);
hold on;

for run = 1:2
    r    = results_coh(run);
    sig  = r.mask;
    cd   = r.coh_diff;
    c    = cols(run,:);

    % Full profile in light colour
    plot(x, cd, '-', 'Color', [c 0.3], 'LineWidth', 1.0, 'HandleVisibility', 'off');

    % Significant segments in solid colour
    seg_plotted = false;
    for i = 1:length(x)-1
        if sig(i) && sig(i+1)
            if ~seg_plotted
                plot(x(i:i+1), cd(i:i+1), '-', 'Color', c, 'LineWidth', 2.5, ...
                    'DisplayName', r.label);
                seg_plotted = true;
            else
                plot(x(i:i+1), cd(i:i+1), '-', 'Color', c, 'LineWidth', 2.5, ...
                    'HandleVisibility', 'off');
            end
        end
    end
    if ~seg_plotted
        % No significant segments — still add to legend
        plot(nan, nan, '-', 'Color', c, 'LineWidth', 2.5, 'DisplayName', r.label);
    end

    % Peak marker
    plot(r.peak_pos, r.coh_diff(r.peak_idx), 'v', ...
        'Color', c, 'MarkerFaceColor', c, 'MarkerSize', 9, ...
        'HandleVisibility', 'off');
end

% Annotate peak positions
for run = 1:2
    r = results_coh(run);
    text(r.peak_pos, r.coh_diff(r.peak_idx), ...
        sprintf('  %.0f mm', r.peak_pos), ...
        'Color', cols(run,:), 'FontSize', 12, 'VerticalAlignment', 'bottom');
end

yline(0, ':k', 'HandleVisibility', 'off');
xlabel('Position along cord (mm, inferior → superior)', 'FontSize', 14);
ylabel('Coherence difference (stat − rest)', 'FontSize', 14);
title(sprintf('Participant 1 — effect of %.1f mm sensor shift on coherence', chosen_mag), ...
    'FontSize', 14, 'FontWeight', 'normal');
legend('Location', 'northwest', 'FontSize', 13, 'Box', 'off');
set(gca, 'FontSize', 14); grid on; box on;

% Print peak shift summary
fprintf('\n── Peak location summary ──\n');
fprintf('  Original:  %.1f mm\n', results_coh(1).peak_pos);
fprintf('  Shifted:   %.1f mm\n', results_coh(2).peak_pos);
fprintf('  Difference: %.1f mm\n', abs(results_coh(2).peak_pos - results_coh(1).peak_pos));

print(gcf, fullfile(fig_dir, ...
    sprintf('coherence_shift_sensitivity_%s.png', chosen_label)), '-dpng', '-r300');
fprintf('Figure saved.\n');

%% =========================================================================
%  LOCAL FUNCTIONS  (copied from main script)
%% =========================================================================
function [coh_diff, cohDiff_perm] = extract_coh_diff(source_perm, nsourcepoints, nPerm)
    cohDiff_perm = zeros(nsourcepoints, nPerm);
    for i = 1:nPerm
        cohDiff_perm(:,i) = source_perm.trialA(i).coh - source_perm.trialB(i).coh;
    end
    coh_diff = source_perm.avgA.coh - source_perm.avgB.coh;
end

function thr = compute_threshold(cohDiff_perm, nsourcepoints)
    maxPerm = max(cohDiff_perm, [], 1);
    thr     = prctile(maxPerm, 95);
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