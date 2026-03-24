%% preproc_spinecoh.m
% Preprocessing pipeline for brain-spine coherence analysis.
% Handles OPM and EMG data for static contraction and rest conditions.
%
% WORKFLOW (two-pass for bad channel detection):
%   Pass 1: Set rmchans = 1 for the FIRST participant run only.
%           This opens an interactive PSD plot — manually select bad
%           channels, then save. Bad channel indices are saved to save_dir.
%   Pass 2: Set rmchans = 0 for all subsequent runs. Bad channels are
%           loaded automatically from the saved file.
%
% EDIT USER CONFIG BELOW THEN RUN.
% =========================================================================

clear all;
close all;

%% =========================================================================
%  USER CONFIG
%  =========================================================================

% --- Toolbox paths --------------------------------------------------------
cfg.spm_path      = 'C:\Users\mspedden\Documents\spm';
cfg.bsc_plot_path = 'C:\Users\mspedden\Documents\brainspineconnectivity\plotting';
cfg.bsc_prep_path = 'C:\Users\mspedden\Documents\brainspineconnectivity\preprocessing';

% --- Data root (BIDS-style: <data_root>/sub-<ID>/ses-001/meg) -------------
cfg.data_root = 'C:\Users\mspedden\Documents';

% --- Output directory -----------------------------------------------------
cfg.save_dir  = 'C:\Users\mspedden\Documents\brainspine_save';

% --- Subjects to process --------------------------------------------------
% To process a single subject, use e.g.: cfg.subs = {'OP00212'};
cfg.subs = {'OP00212','OP00213','OP00215','OP00219', ...
            'OP00220','OP00221','OP00224','OP00225','OP00226'};

% --- Bad channel detection ------------------------------------------------
% rmchans = 1: interactive PSD — select bad channels and save.
%              Only needs to be run once per participant (Pass 1).
% rmchans = 0: load previously saved bad channel file (Pass 2, normal use).
cfg.rmchans = 0;

% --- Processing flags -----------------------------------------------------
cfg.HFCflag   = 1;   % Homogeneous field correction (1 = on)
cfg.HB        = 1;   % Heartbeat removal (1 = on)

% --- Filter parameters ----------------------------------------------------
cfg.fband     = [5 45];     % Bandpass (Hz)
cfg.stpband1  = [48 52];    % Notch 1 (Hz)
cfg.stpband2  = [98 102];   % Notch 2 (Hz)

% --- Epoching / cropping --------------------------------------------------
cfg.EpochTime = 1;     % Epoch length (s)
cfg.DataTime  = 115;   % Duration to crop OPM and EMG to (s)

% --- Channel labels -------------------------------------------------------
cfg.brainchan_labels = {'H1','H6','B3','H8','G3','B6','H7','H4','G4','B4', ...
                        'B2','H2','B7','B8','B5','B1'};
cfg.heart_labels     = {'C1','D3','D2','D4'};
cfg.spine_labels     = {'C2','C3','C4','C7','D8','D5','D6','C5','C8','C6','D7', ...
                        'F2','F7','F5','F8','F3','F1','F4','F6','E2','E4','E8','E5', ...
                        'E1','E7','E6','E3','H5','G5','G7','H3','G1','G2','A7','G6', ...
                        'G8','A6','A3','A2','A8','A5','A1','A4'};
cfg.trigger_label    = 'T5';
cfg.EMGchanname      = 'EXG1';
cfg.chan2replace     = 'Data_Valid2';

% =========================================================================
%  END OF USER CONFIG — do not edit below this line
% =========================================================================

%% Setup
addpath(cfg.spm_path)
spm('defaults','EEG')
addpath(cfg.bsc_plot_path)
addpath(cfg.bsc_prep_path)

if ~exist(cfg.save_dir, 'dir'), mkdir(cfg.save_dir); end

analysis = {'static','rest'};

%% =========================================================================
%  PARTICIPANT LOOP
%% =========================================================================

for subIdx = 1:numel(cfg.subs)

    sub = cfg.subs{subIdx};
    fprintf('\n========================================\n');
    fprintf('  Processing participant %s (%d/%d)\n', sub, subIdx, numel(cfg.subs));
    fprintf('========================================\n');

    % --- Paths for this participant ---------------------------------------
    filetemplate = fullfile(cfg.data_root, ['sub-' sub], 'ses-001', 'meg');
    posfile      = fullfile(filetemplate, 'static_001_ar_positions.tsv');
    EMGpath      = fullfile(cfg.data_root, ['sub-' sub], 'ses-001', 'emg');

    % OP00220 has a difficult heart signal — use best-segment estimation
    flag = strcmp(sub, 'OP00220');

    allfilenames = {};
    beatlen      = [];

    %% ---- Condition loop ------------------------------------------------
    for cond = 1:length(analysis)

        [runs, ref_labels, heartidx, badchanlabels] = get_metadata(sub, analysis{cond});

        %% ---- Run loop --------------------------------------------------
        for k = 1:numel(runs.opm)

            fprintf('  [%s | %s | run %s (%d/%d)]\n', ...
                sub, analysis{cond}, runs.opm{k}, k, numel(runs.opm));

            % Load OPM data
            S     = [];
            S.data      = fullfile(filetemplate, sprintf('%s_%s_array1.lvm', analysis{cond}, runs.opm{k}));
            S.positions = posfile;
            D = spm_opm_create(S);

            %% PSD + bad channel selection
            opms = D.indchannel(D.sensors('MEG').label);

            S           = [];
            S.D         = D;
            S.plot      = 1;
            S.channels  = D.chanlabels(opms);
            S.triallength = 2000;
            S.wind      = @hanning;

            badchan_file = fullfile(cfg.save_dir, sprintf('%s_badchans', sub));

            if cfg.rmchans && k == 1 && cond == 1
                % Pass 1: interactive bad channel selection
                S.selectbad = 1;
                [~,~,badidx] = spm_opm_psd(S);
                if ~isempty(badchanlabels)
                    badsens = find(contains(D.chanlabels, badchanlabels));
                    badidx  = [badidx badsens];
                end
                save(badchan_file, 'badidx');
                fprintf('    Bad channels saved: %s\n', badchan_file);
            else
                % Pass 2: load previously identified bad channels
                S.selectbad = 0;
                spm_opm_psd(S);
                assert(exist([badchan_file '.mat'], 'file') == 2, ...
                    ['Bad channel file not found for %s.\n' ...
                     'Run first with rmchans=1 to generate it:\n  %s'], ...
                    sub, badchan_file);
                load(badchan_file, 'badidx');
            end

            D = badchannels(D, badidx, 1);

            %% Load EMG
            filename = sprintf('%s_00%s.vhdr', sub, runs.emg{k});
            cfg_ft          = [];
            cfg_ft.dataset  = fullfile(EMGpath, filename);
            EMGdata = ft_preprocessing(cfg_ft);

            EMG_trigs = ft_read_event(fullfile(EMGpath, sprintf('%s_00%s.eeg', sub, runs.emg{k})));
            D_EMG     = spm_eeg_ft2spm(EMGdata, sprintf('EMG_%s', sub));

            %% HFC
            if cfg.HFCflag
                fprintf('    Running HFC\n');
                S        = [];
                S.prefix = 'hfc';
                S.D      = D;
                [D_hfc, ~] = spm_opm_hfc(S);

                S          = [];
                S.D1       = D_hfc;
                S.D2       = D;
                S.plot     = 1;
                S.triallength = 2000;
                S.wind     = @hanning;
                spm_opm_rpsd(S);

                D = D_hfc;
            end

            % Mark ref channels as bad so they are excluded from source analysis
            if ~isempty(ref_labels)
                refidx = find(contains(D.chanlabels, ref_labels));
                D      = badchannels(D, refidx, 1);
            end

            %% Filter
            % High-pass
            S        = [];
            S.D      = D;
            S.freq   = min(cfg.fband);
            S.order  = 5;
            S.band   = 'high';
            S.prefix = sprintf('hi%d', cfg.fband(2));
            D        = spm_eeg_filter(S);
            S.D      = D_EMG;
            D_EMG    = spm_eeg_filter(S);

            % Low-pass
            S        = [];
            S.D      = D;
            S.freq   = max(cfg.fband);
            S.order  = 5;
            S.band   = 'low';
            S.prefix = sprintf('lo%d', cfg.fband(2));
            D        = spm_eeg_filter(S);
            S.D      = D_EMG;
            D_EMG    = spm_eeg_filter(S);

            % Notch 1
            S        = [];
            S.D      = D;
            S.band   = 'stop';
            S.order  = 3;
            S.freq   = cfg.stpband1;
            D        = spm_eeg_filter(S);
            S.D      = D_EMG;
            D_EMG    = spm_eeg_filter(S);

            % Notch 2
            S        = [];
            S.D      = D;
            S.band   = 'stop';
            S.order  = 3;
            S.freq   = cfg.stpband2;
            D        = spm_eeg_filter(S);
            S.D      = D_EMG;
            D_EMG    = spm_eeg_filter(S);

            %% Heartbeat estimation
            heart_labels_1 = labelToIndex(cfg.heart_labels);
            heartind = find(contains(D.chanlabels, heart_labels_1));
            megind   = setdiff(D.indchantype('MEG'), D.badchannels);

            if cfg.HB
                if k == 1 && cond == 1
                    [heartest, beatlen] = est_heartbeat(D, heartind, megind, heartidx, flag);
                    allheart = zeros(length(analysis), numel(runs.opm), length(megind), beatlen);
                else
                    heartest = est_heartbeat(D, heartind, megind, heartidx, flag, beatlen);
                end
                allheart(cond, k, :, :) = heartest;
            end

            %% Downsample OPM to match EMG sample rate
            S             = [];
            S.fsample_new = fsample(D_EMG);
            S.D           = D;
            opm_crop      = spm_eeg_downsample(S);

            % Synthetic gradiometer (if ref channels present)
            if ~isempty(ref_labels)
                refidx   = find(contains(opm_crop.chanlabels, ref_labels));
                opm_crop = chantype(opm_crop, refidx, 'REF');
                S        = [];
                S.D      = opm_crop;
                S.confounds = {'REF'};
                opm_crop = spm_opm_synth_gradiometer(S);
            end

            %% Trigger detection and crop
            trigger_idx = find(contains(D.chanlabels, cfg.trigger_label));
            eventCont   = D(trigger_idx, :);
            sample      = 1:length(D.time);
            eventDisc   = (eventCont > 0.9);
            tmp         = eventDisc - [0, eventDisc(1:end-1)];
            events      = round(sample(tmp == 1)');

            % Adjust trigger index for downsampled OPM
            event = round(events / (D.fsample / opm_crop.fsample));

            % Crop OPM
            S        = [];
            S.D      = opm_crop;
            S.timewin = [opm_crop.time(event(1))*1000, ...
                         (opm_crop.time(event(1)) + cfg.DataTime)*1000];
            opm_crop = spm_eeg_crop(S);

            % Crop EMG (second trigger is data start in this system)
            EMG_start = D_EMG.time(EMG_trigs(2).sample);
            S         = [];
            S.D       = D_EMG;
            S.timewin = [EMG_start*1000, (EMG_start + cfg.DataTime)*1000];
            emg_crop  = spm_eeg_crop(S);

            %% Add EMG channel to OPM object
            replaceidx              = opm_crop.indchannel(cfg.chan2replace);
            opm_crop(replaceidx,:)  = 0;
            opm_crop(replaceidx,:)  = emg_crop(1,:);
            opm_crop = chanlabels(opm_crop, replaceidx, cfg.EMGchanname);
            opm_crop = units(opm_crop, replaceidx, 'uV');
            opm_crop = chantype(opm_crop, replaceidx, 'EMG');
            opm_crop.save;

            %% Epoch into fixed-length trials
            epochsamples = cfg.EpochTime * opm_crop.fsample;
            nsamples     = size(opm_crop(:,:), 2);
            Nepochs      = floor(nsamples / epochsamples);

            trl = zeros(Nepochs, 3);
            for f = 1:Nepochs
                trl(f,:) = [(f-1)*epochsamples+1, f*epochsamples, 0];
            end

            condition_labels      = cell(Nepochs, 1);
            condition_labels(:)   = {analysis{cond}};

            S                  = [];
            S.D                = opm_crop;
            S.trl              = trl;
            S.bc               = 0;
            S.conditionlabels  = condition_labels;
            S.prefix           = sprintf('e%dms', round(cfg.EpochTime*1000));
            Dep = spm_eeg_epochs(S);

            %% Outlier trial rejection
            S   = [];
            S.D = Dep;
            [Dep, retain] = spm_opm_removeOutlierTrials(S);

            savename = sprintf('retainedtrials%srun%s_%s', sub, runs.opm{k}, analysis{cond});
            save(fullfile(cfg.save_dir, savename), 'retain');

            allfilenames = strvcat(allfilenames, Dep.fname);

            close all

        end % run loop
    end % condition loop

    %% ---- Merge all runs/conditions for this participant ----------------
    cd(opm_crop.path);
    S        = [];
    S.D      = allfilenames;
    S.prefix = 'merged';
    DAll     = spm_eeg_merge(S);

    %% ---- Heartbeat removal ---------------------------------------------
    if cfg.HB
        avheart = squeeze(mean(allheart, [1 2]));   % average over conditions and runs
        figure; plot(avheart');
        title(sprintf('%s — average heartbeat template', sub), 'Interpreter','none');

        BALANCE      = 1;
        DAll = remove_heartbeat(DAll, avheart, megind, BALANCE);
    end

    fprintf('  Participant %s complete. Merged file: %s\n', sub, DAll.fname);

end % participant loop

fprintf('\n=== PREPROCESSING COMPLETE ===\n');