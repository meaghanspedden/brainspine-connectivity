function XY_coh= cohxy(X,Y,fdat)


[taptr, nfreqs]=size(Y);

%put into fieldtrip structure
Coh_ft=fdat; %copy structure
Coh_ft.label={'Spinal cord'; 'Y'};
Coh_ft.fourierspctrm=zeros(taptr,2,nfreqs);
Coh_ft.fourierspctrm(:,1,:)=X;
Coh_ft.fourierspctrm(:,2,:)=Y;

cfg            = [];
cfg.method     = 'coh';
cfg.channelcmb={'Spinal cord' 'Y'};
XY_coh      = ft_connectivityanalysis(cfg, Coh_ft);
    