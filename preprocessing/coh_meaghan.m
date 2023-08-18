function [cont_coh]=coh_meaghan(D,channels,seed,freqroi)

dat=spm2fieldtrip(D);

% FFT
cfg            = [];
cfg.output     = 'fourier';
cfg.method     = 'mtmfft';
cfg.foilim     = freqroi;
cfg.tapsmofrq  = 5;
cfg.keeptrials = 'yes';


fdat   = ft_freqanalysis(cfg,dat);



tmpcmb=repmat(seed, length(channels),1);
brainchans=channels';


combs=[brainchans tmpcmb];
cfg            = [];
cfg.method     = 'coh';
cfg.channelcmb = combs;

cont_coh      = ft_connectivityanalysis(cfg, fdat);


% chanind=find(contains(cont_coh.labelcmb,'DM-Y'));
% figure
% plot(cont_coh.freq,cont_coh.cohspctrm(chanind,:))

% rightbrain={'19','1B', 'MW', 'MY','OK','MZ'};
% leftbrain={'35','DI','MU','DQ','OH','DG','A1'};


cfg                  = [];
cfg.parameter        = 'cohspctrm';
cfg.xlim             = freqroi;
%cfg.channel=cont_coh.labelcmb(contains(cont_coh.labelcmb,rightbrain));
cfg.refchannel       = seed;
cfg.interplimits     ='electrodes';
cfg.layout           = 'ordered';
cfg.showlabels       = 'yes';
cfg.interactive='yes';
cfg.outline='yes';

%ft_multiplotER(cfg,cont_coh)
