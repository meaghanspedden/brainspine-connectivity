%% Preprocess brain-spine OPM data for FREQ-NESS connectivity analysis
% Pipeline: HFC -> SSP cardiac -> Bandpass -> Notch -> Downsample -> Synth gradiometer
% July 2025

clear all;
close all;

subs = {'OP00213','OP00215', 'OP00219', ...
    'OP00225', 'OP00221', 'OP00224'};

rmchans = 0; % remove bad channels from PSD - set 1 for first pass

analysis = {'static', 'rest'};
save_dir = 'C:\Users\mspedden\Documents\brainspine_save';

if ~exist(save_dir, 'dir')
    mkdir(save_dir);
end

%% Paths
addpath('C:\Users\mspedden\Documents\spm')
spm('defaults','EEG')
addpath('C:\Users\mspedden\Documents\brainspineconnectivity\plotting')
addpath('C:\Users\mspedden\Documents\brainspineconnectivity\preprocessing')

% SSP functions - update path as needed
addpath('C:\Users\mspedden\Documents\brainspineconnectivity\ssp')

%% Channel labels
brainchan_labels = {'H1', 'H6', 'B3', 'H8', 'G3', 'B6', 'H7', 'H4', 'G4', 'B4', 'B2', 'H2', 'B7', 'B8', 'B5', 'B1'};
heart_labels     = {'C1', 'D3', 'D2', 'D4'};

spine_labels = {'C2', 'C3', 'C4', 'C7', 'D8', 'D5', 'D6', 'C5', 'C8', 'C6', 'D7', ...
    'F2', 'F7', 'F5', 'F8', 'F3', 'F1', 'F4', 'F6', 'E2', 'E4', 'E8', 'E5', ...
    'E1', 'E7', 'E6', 'E3', 'H5', 'G5', 'G7', 'H3', 'G1', 'G2', 'A7', 'G6', ...
    'G8', 'A6', 'A3', 'A2', 'A8', 'A5', 'A1', 'A4'};
trigger_label = 'T5';

%% Processing options
fband    = [5 45];
stpband1 = [48 52];
stpband2 = [98 102];
DataTime = 115; % seconds to crop to

HFCflag = 1; % homogeneous field correction

%% Main loop
for s = 1:length(subs)

    sub = subs{s};
    filetemplate = fullfile('C:\Users\mspedden\Documents', ['sub-' sub], 'ses-001', 'meg');
    posfile      = fullfile(filetemplate, 'static_001_ar_positions.tsv');

    for cond = 1:length(analysis)

        [runs, ref_labels, heartidx, badchanlabels] = get_metadata(sub, analysis{cond});

        for k = 1%:numel(runs.opm)

            %% Load data
            S = [];
            S.data      = fullfile(filetemplate, sprintf('%s_%s_array1.lvm', analysis{cond}, runs.opm{k}));
            S.positions = posfile;
            D = spm_opm_create(S);

            %% PSD and bad channel detection
            opms = D.indchannel(D.sensors('MEG').label);

            S         = [];
            S.D       = D;
            S.plot    = 1;
            S.channels   = D.chanlabels(opms);
            S.triallength = 2000;
            S.wind    = @hanning;

            if rmchans && k == 1 && cond == 1
                S.selectbad = 1;
                [~, ~, badidx] = spm_opm_psd(S);

                if ~isempty(badchanlabels)
                    badsens = find(contains(D.chanlabels, badchanlabels));
                    badidx  = [badidx badsens];
                end

                save(fullfile(save_dir, sprintf('%s_badchans', sub)), 'badidx')
            else
                S.selectbad = 0;
                spm_opm_psd(S);
                load(fullfile(save_dir, sprintf('%s_badchans', sub)), 'badidx')
            end

            D = badchannels(D, badidx, 1);

            %% Prepare layout (for visualisation only, kept for compatibility)
            grad = D.sensors('meg');
            grad.coordsys = 'neuromag';

            cfg            = [];
            cfg.grad       = grad;
            cfg.projection = 'orthographic';
            cfg.viewpoint  = 'inferior';
            lay = ft_prepare_layout(cfg);

            %% HFC
            if HFCflag
                S         = [];
                S.prefix  = 'hfc';
                S.D       = D;
                fprintf('\nRunning HFC...\n');
                [D_hfc, ~] = spm_opm_hfc(S);

                % QC: relative PSD before/after HFC
                S         = [];
                S.D1      = D_hfc;
                S.D2      = D;
                S.plot    = 1;
                S.triallength = 2000;
                S.wind    = @hanning;
                spm_opm_rpsd(S);

                D = D_hfc;
            end

            % Mark ref channels as bad so they are excluded from SSP and 
            % source recon, but keep them in the file for synth gradiometry later
            if ~isempty(ref_labels)
                refidx = find(contains(D.chanlabels, ref_labels));
                D = badchannels(D, refidx, 1);
            end

            %% SSP cardiac removal
            % megind: all good MEG channels (excludes bad channels and ref channels)
            megind_ssp = setdiff(D.indchantype('MEG'), D.badchannels);

            % heartind: indices into D of the cardiac monitor channels
            heartind = D.indchannel(heart_labels);

            fprintf('\nRunning SSP cardiac removal (subject %s, cond %s, run %d)...\n', ...
                sub, analysis{cond}, k);

            % Estimate heartbeat template across all good MEG channels
            % heartidx=1: use first cardiac channel for R-peak detection
            % flag=1: find best 30s window for robust R-peak detection
            [heartep, beatlen] = est_heartbeat(D, heartind, megind_ssp, 1, 1, []);

            % Project cardiac components out of continuous data
            % BALANCE=0: do not update forward model (source recon uses VEs downstream)
            D = remove_heartbeat(D, heartep, megind_ssp, 0);

            % Prefix from remove_heartbeat is 'p' (set in hfo_project_out_comps)
            % D is now the SSP-cleaned file

            %% Bandpass filter
            S        = [];
            S.D      = D;
            S.freq   = min(fband);
            S.order  = 5;
            S.band   = 'high';
            S.prefix = sprintf('hi%d', fband(1));
            D = spm_eeg_filter(S);

            S        = [];
            S.D      = D;
            S.freq   = max(fband);
            S.order  = 5;
            S.band   = 'low';
            S.prefix = sprintf('lo%d', fband(2));
            D = spm_eeg_filter(S);

            %% Notch filters
            S        = [];
            S.D      = D;
            S.band   = 'stop';
            S.order  = 3;
            S.freq   = stpband1;
            D = spm_eeg_filter(S);

            S        = [];
            S.D      = D;
            S.band   = 'stop';
            S.order  = 3;
            S.freq   = stpband2;
            D = spm_eeg_filter(S);

            %% Downsample
            S            = [];
            S.fsample_new = 125;
            S.D          = D;
            opm_crop = spm_eeg_downsample(S);

            %% Synthetic gradiometry (reference channel denoising)
            % Unmark ref channels so spm_opm_synth_gradiometer can use them,
            % then re-mark after
            if ~isempty(ref_labels)
                refidx_crop = find(contains(opm_crop.chanlabels, ref_labels));
                opm_crop    = chantype(opm_crop, refidx_crop, 'REF');

                S           = [];
                S.D         = opm_crop;
                S.confounds = {'REF'};
                opm_crop    = spm_opm_synth_gradiometer(S);
            end

            fprintf('\nPreprocessing complete: subject %s, condition %s, run %d\n', ...
                sub, analysis{cond}, k);

        end % runs
    end % conditions
end % subjects