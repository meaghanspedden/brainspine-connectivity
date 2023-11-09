function [cspect,fdat,trcspect]=xspectrum_meaghan(D,channels,seed,freqroi,seedpermute,normseed,replace)

if nargin<5
    seedpermute=[];
end;
if nargin<6
    normseed=0;
end;
if nargin<7,
    replace=[];
end;

if isempty(seedpermute),
    seedpermute=0;
end;
if isempty(normseed),
    normseed=0;
end;

dat=spm2fieldtrip(D);

seedind=find(contains(dat.label,seed));
if seedpermute,
    warning('Permuting seed chan trials')
    Ntrials=length(dat.trial);
    rtrial=randperm(Ntrials);
    
    for f=1:Ntrials
        dum=dat.trial{f};
        dum2=dat.trial{rtrial(f)};
        trialdata=dum;
        trialdata(seedind,:)=dum2(seedind,:);
        dat.trial{f}=trialdata; 
    end
end % if
    


% FFT
cfg            = [];
cfg.output     = 'fourier';
cfg.method     = 'mtmfft';
cfg.foilim     = freqroi;
cfg.tapsmofrq  = 5;
cfg.keeptrials = 'yes';
cfg.channel=[channels,seed];


fdat   = ft_freqanalysis(cfg,dat);

if ~isempty(replace),
    fprintf('\n Replacing %s with linear mix of channels\n',replace.label)
    ind=find(contains(fdat.label,replace.label));
    for t=1:size(fdat.fourierspctrm,1)
        data=squeeze(fdat.fourierspctrm(t,replace.ind,:));
        ndata=replace.U'*data;
        fdat.fourierspctrm(t,ind,:)=ndata;
    end;
end;

if normseed,
    warning('Normalizing seed chan trials')
    fspect=fdat.fourierspctrm(:,seedind,:);
    absfspect=sqrt(fspect.*conj(fspect));
    fspect=fspect./absfspect;
    fdat.fourierspctrm(:,seedind,:)=fspect;
end;

tmpcmb=repmat(seed, length(channels),1);
brainchans=channels';


combs=[brainchans tmpcmb];
cfg            = [];
cfg.method     = 'csd';
cfg.complex='complex';
cfg.channelcmb = combs;


[cspect,trialdata]      = ft_connectivityanalysis_grb_v2(cfg, fdat); %% GRB hack to get trialdata from this

trcspect=trialdata.crsspctrm(:,1:length(channels),:);
 







