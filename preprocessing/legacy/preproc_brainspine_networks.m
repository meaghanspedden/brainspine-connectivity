%% process data from hand contraction. ~ July 2025


clear all;
close all;

subs = {'OP00213','OP00215', 'OP00219', ...
    'OP00225', 'OP00221', 'OP00224'}; %have brain

rmchans=0; %remove bad channels from psd...for first round of analysis

analysis={'static', 'rest'};
save_dir = 'C:\Users\mspedden\Documents\brainspine_save';

if ~exist(save_dir, 'dir')
    mkdir(save_dir);
end


%% paths
addpath('C:\Users\mspedden\Documents\spm')
spm('defaults','EEG')
addpath('C:\Users\mspedden\Documents\brainspineconnectivity\plotting')
addpath('C:\Users\mspedden\Documents\brainspineconnectivity\preprocessing')


%% labels
brainchan_labels={'H1', 'H6', 'B3', 'H8', 'G3', 'B6', 'H7', 'H4', 'G4', 'B4', 'B2', 'H2', 'B7', 'B8', 'B5', 'B1'};
heart_labels={'C1', 'D3', 'D2', 'D4'};

spine_labels={'C2', 'C3', 'C4', 'C7', 'D8', 'D5', 'D6', 'C5', 'C8', 'C6', 'D7', ...
    'F2', 'F7', 'F5', 'F8', 'F3', 'F1', 'F4', 'F6', 'E2', 'E4', 'E8', 'E5', ...
    'E1', 'E7', 'E6', 'E3', 'H5', 'G5', 'G7', 'H3', 'G1', 'G2', 'A7', 'G6', ...
    'G8', 'A6', 'A3', 'A2', 'A8', 'A5', 'A1', 'A4'};
trigger_label='T5';



%% Processing Options
fband=[5 45]; %% filter band
stpband1=[48 52];
stpband2=[98 102];
DataTime=115; %length to crop opm and emg to, seconds

% flags to turn these filters on/off
HFCflag=1; %% homogenous field correction (0 or 1)


for s=1:length(subs)

    sub=subs{s};

    filetemplate = fullfile('C:\Users\mspedden\Documents', ['sub-' sub], 'ses-001', 'meg');

    posfile=fullfile(filetemplate, 'static_001_ar_positions.tsv');
    %% preproc loop

    allfilenames={};

    for cond=1:length(analysis)

        [runs, ref_labels, heartidx, badchanlabels] = get_metadata(sub,analysis{cond}); %only runs differ here woul be better to output across condition


        for k=1%:numel(runs.opm)

            S=[];
            S.data = fullfile(filetemplate, sprintf('%s_%s_array1.lvm', analysis{cond},runs.opm{k}));
            S.positions = posfile;
            D = spm_opm_create(S);

            %% psd
            opms=D.indchannel(D.sensors('MEG').label);

            S = [];
            S.D = D;
            S.plot = 1;
            S.channels = D.chanlabels(opms);
            S.triallength = 2000;
            S.wind = @hanning;

            if rmchans && k==1 && cond==1
                S.selectbad=1;
                [~,~,badidx] = spm_opm_psd(S);


                if~isempty(badchanlabels)
                    badsens=find(contains(D.chanlabels,badchanlabels));
                    badidx=[badidx badsens];
                end

                save(fullfile(save_dir,sprintf('%s_badchans',sub)), 'badidx')

            else
                S.selectbad=0;
                spm_opm_psd(S);
                load(fullfile(save_dir,sprintf('%s_badchans',sub)), 'badidx') %same bad channels across runs
            end

            D=badchannels(D, badidx,1);

            grad = D.sensors('meg');
            grad.coordsys = 'neuromag';

            cfg = [];
            cfg.grad = grad;
            cfg.projection = 'orthographic';
            cfg.viewpoint = 'inferior';

            lay = ft_prepare_layout(cfg);

            if HFCflag
                S=[];
                S.prefix='hfc';
                S.D=D;
                %S.L=2;
                fprintf('\n Running HFC');
                [D_hfc, ~]=spm_opm_hfc(S);

                S = [];
                S.D1 =D_hfc;
                S.D2=D;
                S.plot = 1;
                %S.channels = D.chanlabels(setdiff(opms,badidx));
                S.triallength = 2000;
                S.wind = @hanning;
                [shield,~] = spm_opm_rpsd(S);

                D=D_hfc;
            end

            if ~isempty(ref_labels) %ref channels included in hfc then set to bad so as not to include in source analysis
                refidx=find(contains(D.chanlabels,ref_labels));
                D=badchannels(D, refidx,1);

            end
            %%  Filter

            S=[];
            S.D=D;
            S.freq=min(fband);
            S.order=5;
            S.band='high';
            S.prefix=sprintf('hi%d',fband(2));
            D=spm_eeg_filter(S);


            S=[];
            S.D=D;
            S.freq=max(fband);
            S.order=5;
            S.band='low';
            S.prefix=sprintf('lo%d',fband(2));
            D=spm_eeg_filter(S);


            S=[];
            S.D=D;
            S.band='stop';
            S.order=3;
            S.freq=stpband1;
            S.D=D;
            D=spm_eeg_filter(S);

            S=[];
            S.D=D;
            S.band='stop';
            S.order=3;
            S.freq=stpband2;
            D=spm_eeg_filter(S);


            megind=setdiff(D.indchantype('MEG'), D.badchannels);


            %% downsample OPM data

            S=[];
            S.fsample_new=125;
            S.D=D;
            opm_crop=spm_eeg_downsample(S);

            if ~isempty(ref_labels)
                refidx=find(contains(opm_crop.chanlabels,ref_labels));
                opm_crop=chantype(opm_crop,refidx,'REF');

                S=[];
                S.D=opm_crop;
                S.confounds={'REF'};

                opm_crop = spm_opm_synth_gradiometer(S);
            end


            %% currently this cant handle hfc previously
            addpath(genpath('C:\Users\mspedden\Documents\opm_ica'))

            S           = [];
            S.D         = opm_crop;        % meeg object
            S.layout    = lay;     
            S.ncomps    = 50;       % number of components
            S.infer     = true;     % run ICA
            S.review    = true;    % don't do component review

            iD = spm_opm_ica(S);

            %visual inspection of time series
            %         ftdat=spm2fieldtrip(opm_crop);
            %         ftdat=r2mfield(ftdat, 'hdr');
            %         badidx=D.badchannels;
            %         badlabs=ftdat.label(badidx);
            %         labels = ftdat.label;
            %         isGood = ~contains(labels, badlabs);
            %         hasXYZ = contains(labels, 'X') | contains(labels, 'Y') | contains(labels, 'Z');
            %         goodlabs = labels(isGood & hasXYZ);
            %         cfg=[];
            %         cfg.channel=goodlabs;
            %         cfg.allowoverlap='yes';
            %         ft_databrowser(cfg,ftdat)

            %close all
        end % for runs
    end %loop through conditions
end %subjects



