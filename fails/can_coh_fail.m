%% CVA

D=spm_eeg_load('D:\MSST001\sub-OP00189\ses-001\meg\mergeoe1000mse1000msdfflo45hi45sub-OP00189_ses-001_task-static_run-001_meg.mat');

%brainchan_labels
%ref_labels
%sc_labels

emg_label=all_labels(EMGidx);

goodchanidx=setdiff(1:size(D,1),badchannels(D));

ftdat=spm2fieldtrip(D);
ftdat=rmfield(ftdat,'hdr');

cfg=[];
cfg.channel=D.chanlabels(goodchanidx);
ft_selectdata(cfg,ftdat);


cfg            = [];
cfg.output     = 'powandcsd';
cfg.method     = 'mtmfft';
cfg.foilim     = [20 20];
cfg.tapsmofrq  = 1;
cfg.keeptrials = 'no';
cfg.channel=brainchan_labels(contains(brainchan_labels, '-Z'));

freq_dat   = ft_freqanalysis(cfg, ftdat);

%% reformat output
% Extract the cross-spectrum and channel pairs from the FieldTrip output
crsspctrm = freq_dat.crsspctrm;  % Size: 561x1 (cross-spectra for channel pairs)
labelcmb = freq_dat.labelcmb;    % Size: 561x2 (pairs of channels)

% Initialize the auto-covariance matrix (Cx) to be M x M (34 x 34 in this case)
M = length(freq_dat.label);  % Number of channels (34)
Cx = zeros(M, M);  % Initialize auto-covariance matrix

% Loop through the channel pairs and populate the cross-spectra
for pair_idx = 1:size(labelcmb, 1)
    % Find the indices of the channels in the pair
    ch1_idx = find(strcmp(freq_dat.label, labelcmb{pair_idx, 1}));
    ch2_idx = find(strcmp(freq_dat.label, labelcmb{pair_idx, 2}));
    
    % Assign the cross-spectrum value to the corresponding place in the matrix
    Cx(ch1_idx, ch2_idx) = crsspctrm(pair_idx);
    Cx(ch2_idx, ch1_idx) = crsspctrm(pair_idx);  % Matrix is symmetric
end

% Cx is n sc chans x n sc chans, complex (for y only)
%%

cfg            = [];
cfg.output     = 'fourier';
cfg.method     = 'mtmfft';
cfg.foilim     = [20 20];
cfg.tapsmofrq  = 1;
cfg.keeptrials = 'yes';
cfg.channel=ftdat.label(contains(ftdat.label,EMGchan));

freq_emg   = ft_freqanalysis(cfg, ftdat);

Cy=mean(freq_emg.fourierspctrm);


%%

seed=emg_label;
megchans=brainchan_labels(contains(brainchan_labels, '-Z'))';
tmpcmb=repmat(seed, length(megchans),1);
combs=[megchans tmpcmb];
cfg            = [];
cfg.output     = 'powandcsd';
cfg.method     = 'mtmfft';
cfg.foilim     = [20 20];
cfg.tapsmofrq  = 1;
cfg.keeptrials = 'no';
cfg.channelcmb=combs;

freq_both   = ft_freqanalysis(cfg, ftdat);

Cxy=freq_both.crsspctrm;

[coh_complex, wa, wb, ta, tb] = canonical_coherence(Cx, Cy, Cxy);

cancoh=coh;
cancoh.powspctrm=ta;
cancoh.label=megchans;
cancoh.freq=20;

fig=figure;
fig.Position = [100, 200, 600, 400];
cfg                  = [];
cfg.parameter        = 'powspctrm';
cfg.xlim             = [20 20];
cfg.zlim             ='zeromax';
cfg.interplimits     ='electrodes';
cfg.layout           = lay;

ft_topoplotER(cfg,cancoh)
colorbar

