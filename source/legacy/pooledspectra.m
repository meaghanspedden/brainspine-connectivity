%% pooled_spectra.m
% Pooled NeuroSpec population analysis across subjects.
%
% Computes pooled spectra and coherence for three pairs:
%   Brain VE <-> EMG
%   Brain VE <-> Spine VE
%   Spine VE  <-> EMG
%
% Uses pool_sc / pool_sc_out (NeuroSpec 2.0 population routines).
% Each pair accumulates spectral coefficients (sc) across subjects,
% then pool_sc_out converts the result to plottable f/t/cl structs.
%
% Requires step5 VE files to already exist in save_dir.
%
% Outputs (written to fig_dir):
%   pooled_spectra_brainEMG.fig / .png
%   pooled_spectra_brainSpine.fig / .png
%   pooled_spectra_spineEMG.fig / .png
%   pooled_spectra_all3pairs.fig / .png

clear all; close all; clc;

%% =========================================================================
%  USER CONFIG
%% =========================================================================
fieldtrip_path  = 'C:\Users\mspedden\Documents\fieldtrip';
spm_path        = 'C:\Users\mspedden\Documents\spm';
bsc_path        = 'C:\Users\mspedden\Documents\brainspineconnectivity\source';
neurospec_path  = 'C:\Users\mspedden\Documents\neurospec20';

data_root       = 'C:\spinecoh_data';
save_dir        = 'C:\Users\mspedden\Documents\brainspine_savetest';

brain_ve_suffix  = '_brain_pct10';
spine_ve_pattern = 'VE_spine_prevalence_sub%s_forspectra_BS.mat';
spine_ve_varname = 'VE';

subs = {'OP00212','OP00213','OP00215','OP00219', ...
        'OP00225','OP00221','OP00224'};

seg_pwr = 11;          % segment length = 2^seg_pwr samples
fband   = [10 35];     % beta band for shading
opt_str = 'M4';        % multitaper option

saveFigs = 1;

%% =========================================================================
%  SETUP
%% =========================================================================
addpath(bsc_path); addpath(spm_path);
spm('defaults','EEG');
addpath(fieldtrip_path);
ft_defaults;
addpath(genpath(neurospec_path));

fig_dir = fullfile(save_dir, 'figures', 'pooled_spectra');
if ~exist(fig_dir,'dir'), mkdir(fig_dir); end

nSubs   = numel(subs);
seg_len = 2^seg_pwr;

% Accumulators for pool_sc — initialise as empty (first call creates them)
plf_BE = []; plv_BE = [];   % Brain-EMG
plf_BS = []; plv_BS = [];   % Brain-Spine
plf_SE = []; plv_SE = [];   % Spine-EMG

fprintf('\n=== POOLED SPECTRA ===\n');
fprintf('  Subjects: %d\n', nSubs);
fprintf('  seg_pwr: %d  (seg_len=%d)\n', seg_pwr, seg_len);
fprintf('======================\n\n');

%% =========================================================================
%  SUBJECT LOOP
%% =========================================================================
for ss = 1:nSubs
    sub = subs{ss};
    fprintf('  Subject %s (%d/%d)...\n', sub, ss, nSubs);

    %% Load raw data
    run = '001'; if strcmp(sub,'OP00224'), run = '002'; end
    datafile = fullfile(data_root, ['sub-' sub], 'ses-001', 'meg', ...
        sprintf('pmergedoe1000mspddfflo45hi45hfcstatic_%s_array1.mat', run));
    D     = spm_eeg_load(datafile);
    ftdat = spm2fieldtrip(D);

    %% Load brain VE
    bVE_file = fullfile(save_dir, ...
        sprintf('sub%s_VE_brain_M1%s.mat', sub, brain_ve_suffix));
    bVE = load(bVE_file, 'VE_brain');
    bVE = bVE.VE_brain;

    %% Load spine VE
    sVE_file = fullfile(save_dir, sprintf(spine_ve_pattern, sub));
    sVE_data = load(sVE_file);
    sVE      = sVE_data.(spine_ve_varname);

    %% Rectified EMG
    cfg_ft = []; cfg_ft.channel = 'EXG1';
    EMG = ft_selectdata(cfg_ft, ftdat);
    cfg_ft = []; cfg_ft.rectify = 'yes';
    EMG = ft_preprocessing(cfg_ft, EMG);

    %% Static contraction trials only
    statidx = find(ftdat.trialinfo == 1);
    cfg_ft  = []; cfg_ft.trials = statidx;
    statB   = ft_selectdata(cfg_ft, bVE);
    statS   = ft_selectdata(cfg_ft, sVE);
    statEMG = ft_selectdata(cfg_ft, EMG);

    %% Concatenate into continuous signals
    Bcont   = [statB.trial{:}];
    Scont   = [statS.trial{:}];
    EMGcont = abs([statEMG.trial{:}]);   % rectified
    samp_rate = ftdat.hdr.Fs;

    %% NeuroSpec — compute spectra + spectral coefficients (4th output = sc)
    % sp2a2_R2_mt returns: [f, t, cl, sc]
    %   f  = spectral estimates matrix
    %   t  = cumulant density matrix
    %   cl = confidence limit structure
    %   sc = spectral coefficients (needed for pool_sc)
    [~, ~, cl_BE, sc_BE] = sp2a2_R2_mt(Bcont',   EMGcont', samp_rate, seg_pwr, opt_str);
    [~, ~, cl_BS, sc_BS] = sp2a2_R2_mt(Bcont',   Scont',   samp_rate, seg_pwr, opt_str);
    [~, ~, cl_SE, sc_SE] = sp2a2_R2_mt(Scont',   EMGcont', samp_rate, seg_pwr, opt_str);

    %% Accumulate into pooled analysis
    % First subject: initialise.  Subsequent: update.
    if isempty(plf_BE)
        [plf_BE, plv_BE] = pool_scf(sc_BE, cl_BE);
        [plf_BS, plv_BS] = pool_scf(sc_BS, cl_BS);
        [plf_SE, plv_SE] = pool_scf(sc_SE, cl_SE);
    else
        [plf_BE, plv_BE] = pool_scf(sc_BE, cl_BE, plf_BE, plv_BE);
        [plf_BS, plv_BS] = pool_scf(sc_BS, cl_BS, plf_BS, plv_BS);
        [plf_SE, plv_SE] = pool_scf(sc_SE, cl_SE, plf_SE, plv_SE);
    end
end

fprintf('\n  All subjects accumulated.  Converting to plottable form...\n');

%% =========================================================================
%  CONVERT POOLED RESULTS TO PLOTTABLE FORM
%% =========================================================================
% pool_sc_out returns f, t, cl in same format as sp2a2_R2_mt core output
[f_BE, t_BE, cl_BE_pool] = pool_scf_out(plf_BE, plv_BE);
[f_BS, t_BS, cl_BS_pool] = pool_scf_out(plf_BS, plv_BS);
[f_SE, t_SE, cl_SE_pool] = pool_scf_out(plf_SE, plv_SE);

%% =========================================================================
%  PLOTS — individual pair figures
%% =========================================================================
col_BE = [0.2 0.4 0.8];
col_BS = [0.8 0.2 0.4];
col_SE = [0.2 0.7 0.3];

pairs = {
    f_BE, t_BE, cl_BE_pool, col_BE, 'Brain VE – EMG',       'brainEMG';
    f_BS, t_BS, cl_BS_pool, col_BS, 'Brain VE – Spine VE',  'brainSpine';
    f_SE, t_SE, cl_SE_pool, col_SE, 'Spine VE – EMG',       'spineEMG';
};

for pp = 1:size(pairs,1)
    f_p  = pairs{pp,1};
    t_p  = pairs{pp,2};
    cl_p = pairs{pp,3};
    col  = pairs{pp,4};
    ttl  = pairs{pp,5};
    tag  = pairs{pp,6};

    hfig = figure('Color','w','Position',[100 100 900 700]);

    % --- Row 1: input / output log spectra ---
    subplot(3,2,1);
    plot_pooled_spectrum(f_p(:,1), f_p(:,2), cl_p.f_c95, fband, col, ...
        'Input log-spectrum', 'Log power');

    subplot(3,2,2);
    plot_pooled_spectrum(f_p(:,1), f_p(:,3), cl_p.f_c95, fband, col, ...
        'Output log-spectrum', 'Log power');

    % --- Row 2: coherence ---
    subplot(3,2,[3 4]);
    fax = f_p(:,1);
    coh = f_p(:,4);
    conf95 = cl_p.ch_c95;
    hold on;
    fill_band(fax, fband, [0.93 0.93 0.93]);
    yline(conf95, '--k', 'LineWidth', 1.2, 'HandleVisibility','off');
    plot(fax, coh, '-', 'Color', col, 'LineWidth', 2);
    xlabel('Frequency (Hz)', 'FontSize', 12);
    ylabel('Coherence', 'FontSize', 12);
    title(sprintf('Pooled coherence  (n=%d)  |  95%% CL = %.4f', nSubs, conf95), ...
        'FontSize', 12);
    xlim([0 60]); set(gca,'FontSize',11); box off;

    % --- Row 3: cumulant density ---
    subplot(3,2,[5 6]);
    tms  = t_p(:,1);   % lag axis (ms)
    cum  = t_p(:,2);   % cumulant
    conf_t = cl_p.q_c95;
    hold on;
    yline( conf_t, '--k', 'LineWidth', 1.2, 'HandleVisibility','off');
    yline(-conf_t, '--k', 'LineWidth', 1.2, 'HandleVisibility','off');
    yline(0, '-k', 'LineWidth', 0.8, 'HandleVisibility','off');
    plot(tms, cum, '-', 'Color', col, 'LineWidth', 1.8);
    xlabel('Lag (ms)', 'FontSize', 12);
    ylabel('Cumulant density', 'FontSize', 12);
    title('Pooled cumulant density', 'FontSize', 12);
    xlim([-100 100]); set(gca,'FontSize',11); box off;

    sgtitle(sprintf('Pooled spectra — %s  (n=%d subjects)', ttl, nSubs), ...
        'FontSize', 14, 'FontWeight','bold', 'Interpreter','none');

    if saveFigs
        savefig(hfig, fullfile(fig_dir, sprintf('pooled_spectra_%s.fig', tag)));
        saveas(hfig,  fullfile(fig_dir, sprintf('pooled_spectra_%s.png', tag)));
    end
end

%% =========================================================================
%  SUMMARY FIGURE — coherence only, all 3 pairs on one axes
%% =========================================================================
hfig_all = figure('Color','w','Position',[100 100 700 400]);
hold on;
fill_band(f_BE(:,1), fband, [0.93 0.93 0.93]);

plot(f_BE(:,1), f_BE(:,4), '-', 'Color', col_BE, 'LineWidth', 2, 'DisplayName','Brain–EMG');
plot(f_BS(:,1), f_BS(:,4), '-', 'Color', col_BS, 'LineWidth', 2, 'DisplayName','Brain–Spine');
plot(f_SE(:,1), f_SE(:,4), '-', 'Color', col_SE, 'LineWidth', 2, 'DisplayName','Spine–EMG');

% Confidence limits (use Brain-EMG cl as representative; all same n)
yline(cl_BE_pool.ch_c95, '--k', 'LineWidth', 1.2, 'HandleVisibility','off');

xlabel('Frequency (Hz)', 'FontSize', 13);
ylabel('Coherence', 'FontSize', 13);
title(sprintf('Pooled coherence — all pairs  (n=%d)', nSubs), ...
    'FontSize', 14, 'Interpreter','none');
legend('Location','northeast', 'Box','off', 'FontSize', 12);
xlim([0 60]); set(gca,'FontSize',12); box off;

if saveFigs
    savefig(hfig_all, fullfile(fig_dir, 'pooled_spectra_all3pairs.fig'));
    saveas(hfig_all,  fullfile(fig_dir, 'pooled_spectra_all3pairs.png'));
end

fprintf('\n=== POOLED SPECTRA DONE ===\n');
fprintf('  Figures saved to: %s\n\n', fig_dir);

%% =========================================================================
%  LOCAL HELPERS
%% =========================================================================

function plot_pooled_spectrum(fax, logpow, conf95, fband, col, ttl, ylbl)
    hold on;
    fill_band(fax, fband, [0.93 0.93 0.93]);
    plot(fax, logpow, '-', 'Color', col, 'LineWidth', 1.8);
    xlabel('Frequency (Hz)', 'FontSize', 11);
    ylabel(ylbl, 'FontSize', 11);
    title(ttl, 'FontSize', 11);
    xlim([0 60]); set(gca,'FontSize',10); box off;
end

function fill_band(fax, fband, col)
    yl = ylim;
    fill([fband(1) fband(2) fband(2) fband(1)], ...
         [yl(1) yl(1) yl(2) yl(2)], col, ...
         'FaceAlpha', 0.5, 'EdgeColor','none', 'HandleVisibility','off');
end