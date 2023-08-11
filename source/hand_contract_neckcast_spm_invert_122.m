%% SPM inversion - new neck cast - hand contraction July 2023

clear all;
close all;

mydir='D:\MSST001';


%filenames for concatenated runs for right hand and left hand
filenames=strvcat([mydir '\sub-OP00122\ses-001\meg\pprhandoe1000msfdfflo45hi5ds_sub-OP00122_ses-001_task-static_right_run-002_meg.mat'],...
        [mydir '\sub-OP00122\ses-001\meg\pplhandoe1000msfdfflo45hi5ds_sub-OP00122_ses-001_task-static_left_run-001_meg.mat']);

posfile=[mydir '\sub-OP00122\ses-001\meg\ds_sub-OP00122_ses-001_positions.tsv'];
backshape='D:\OP00122_experiment\cast space\OP0015_seated.stl'; %stl in same space as sensors
cyl='D:\OP00122_experiment\cast space\cylinder_good_space.stl'; %cylinder for source space

dpath=[mydir '\sub-OP00122\ses-001\meg']; %data path
savepath=[mydir '\Coh_results00122']; %save path

if ~exist(savepath,'dir')
    mkdir(savepath)
end

cd(dpath)

brainchan_labels={'19','DG', 'OH', 'A1','1B', 'A9','JS','A6','DJ','MY','DS','OK','MI','17'};
badchans={'ML-X','ML-Y','ML-Z','K4-Z','K4-X'}; %% these are also marked in D object

%%
addpath D:\torso_tools
addpath D:\spm12
spm('defaults','EEG');
addpath(genpath('D:\brainspineconnectivity'))

%% analysis options
SHUFFLE=1;
rng(125) %for shuffling
allcanfilenames=[];

whatstr='emg+abs';
%whatstr='brainopt+abs';
%whatstr='orthbrain+brainopt';


FIXORIENT=[]; 
REMOVELIN=1; %this seems to be hard coded now
RMORTHBRAIN=contains(whatstr,'brain'); %% don't remove brain if dealing with cord-muscle
FIXMODES=0;
CLONETRIALS=1; % use all trial data (1) rather than average (0)

freqroi=[5 35];

%invtype='EBBr95';
invtype='IID';
%invtype='EBB';

IMAGECROSS=0;
COMB3AXES=1;


%colors for lines in plots
cols=colormap(brewermap([],"Dark2"));
col1=cols(1,:); col2=cols(4,:); %for non linear
col3=cols(3,:); col4=cols(6,:); %for real imag and coh
col5=cols(5,:);



for cnd=1%:size(filenames,1),

    DAll=spm_eeg_load(filenames(cnd,:));
    [a1,b1,c1]=fileparts(DAll.fullfile);
    pstr=b1(1:10); %%get upper/lower l/r info from filename

    cname=[a1 filesep 'clone_b1' b1 c1 ];

    D=DAll;

    %% identify channels on spine (the ones that have positions and ori)

    grad=D.sensors('MEG');
    badgrad=badchannels(D);
    [a1,b1,spineind]=intersect(grad.label,D.chanlabels);
    spineind=setdiff(spineind,D.badchannels); %% remove bad channels


    %% get list of indices of channels on cortex
    brainind=[];

    for g=1:length(brainchan_labels)
        brainind=[brainind find(contains(D.chanlabels,deblank(brainchan_labels(g))))];
    end

    emgind=D.indchantype('EMG');


    %% now compute coherence and cross spectra between brain, emg, and sensors on cord
    cd  D:\fieldtrip-master\fieldtrip-master\private %fix this?

    seedperm=0;
    cohbrain=coh_meaghan(D,D.chanlabels(brainind),D.chanlabels(emgind(1)),[min(freqroi) max(freqroi)]);
    [cspect_brain,fdat_brain,trialcspect_brain]=xspectrum_meaghan(D,D.chanlabels(brainind),D.chanlabels(emgind(1)),[min(freqroi) max(freqroi)],seedperm);
    [cspect,fdat,trialcspect]=xspectrum_meaghan(D,D.chanlabels(spineind),D.chanlabels(emgind(1)),[min(freqroi) max(freqroi)],seedperm);

    fbemgind=find(contains(fdat_brain.label,D.chanlabels(emgind(1))));
    fbrainind=setdiff(1:length(fdat_brain.label),fbemgind);

    %     %% plot some metrics of brain-emg interaction
    %     nchansBrain=length(cspect_brain.labelcmb); %it was getting confused because csd mat is square
    %
    %     figure; subplot(4,1,1);
    %     plot(cspect_brain.freq,real(cspect_brain.crsspctrm(1:nchansBrain,:))');
    %     title('real cross')
    %     subplot(4,1,2);
    %     plot(cspect_brain.freq,imag(cspect_brain.crsspctrm(1:nchansBrain,:))');
    %     title('imag cross')
    %     subplot(4,1,3)
    %     plot(cspect_brain.freq,abs(cspect_brain.crsspctrm(1:nchansBrain,:))');
    %     title('abs cross')
    %     subplot(4,1,4)
    %     plot(cohbrain.freq,cohbrain.cohspctrm(1:nchansBrain,:)');
    %     title('coh')
    %% get dominant spatial component that explains the brain-emg cross spectrum
    r1=cspect_brain.crsspctrm;
    r1=r1-mean(r1);

    [Ur,S,~]=svd(real(r1)*real(r1)'); %% could be an issue if peak coherence is imag.
    varexp=cumsum(diag(S))./sum(diag(S));
    Nrcomp=min(find(varexp>0.99));
    %Ur=Ur(:,1:Nrcomp);
    [Ui,Si,Vi]=svd(imag(r1)*imag(r1)');

    Uinv=pinv(Ur); %% the spatial mixture of brain chans orthogonal to the ideal mixture for EMG-brain cross-spectrum
    %CHECKMIX=0;
    %if CHECKMIX
    replace.U=Ui(:,1); %% Ur or Uinv
    replace.ind=fbrainind;
    replace.label='G2-19-Y';

    %calculate coherence with ideal mixture-EMG
    [~, fdat_brainmix,~]=xspectrum_meaghan(D,D.chanlabels(brainind),D.chanlabels(emgind(1)),[min(freqroi) max(freqroi)],seedperm,[],replace);

    combs={replace.label D.chanlabels{emgind(1)}};
    cfg            = [];
    cfg.method     = 'coh';
    cfg.channelcmb = combs;

    brainmix_emg_coh      = ft_connectivityanalysis(cfg, fdat_brainmix);

    %% replace one of the chans with optimcal linear mixture and check it is larger/smaller
   
    [cspect_brainmix,fdat_brainmix,trialcspect_brainmix]=xspectrum_meaghan(D,D.chanlabels(brainind),D.chanlabels(emgind(1)),[min(freqroi) max(freqroi)],seedperm,[],replace);
    %             figure;
    %             plot(cspect_brain.freq,abs(cspect_brain.crsspctrm),'r');
    %             hold on;
    %             plot(cspect_brain.freq,abs(cspect_brainmix.crsspctrm),'go');
    ind=find(contains(cspect_brainmix.labelcmb,replace.label));
    %             plot(cspect_brain.freq,abs(cspect_brainmix.crsspctrm),'go');
    %             plot(cspect_brain.freq,abs(cspect_brainmix.crsspctrm(ind,:)),'k*');
    %             title('abs cross')
    % end % checkmix

    %% get indices of channels within the fdat structure

    fspemgind=find(contains(fdat.label,D.chanlabels(emgind(1))));
    spind=setdiff(1:length(fdat.label),fspemgind);
    bind=setdiff(1:length(fdat_brain.label),fbemgind);

    femgind=find(contains(fdat.label,D.chanlabels(emgind(1))));
    fuseind=setdiff(1:length(spineind),femgind); %freqs to use

    %% average over taps to get 1 spectrum per trial

    trialfspect=zeros(D.ntrials,length(spineind),size(fdat.fourierspctrm,3));

    for f=1:D.ntrials
        nt=fdat.cumtapcnt(f);
        trialind=(f-1)*nt+1:f*nt;
        trialfspect(f,:,:)=squeeze(mean(fdat.fourierspctrm(trialind,fuseind,:)));
    end


    if ~IMAGECROSS
        warning('just imaging power not cross spectrum');
        trialcspect=trialfspect;
    end

    Nclonesamples=length(cspect.freq)*2; % x2 for real and imag part
    cind=[];
    cind(1,:)=1:(Nclonesamples/2);
    cind(2,:)=((Nclonesamples/2)+1):Nclonesamples;

    imagefreq=cspect.freq;

    if CLONETRIALS
        Dclone=clone(DAll,cname,[size(DAll,1) Nclonesamples size(DAll,3)]);
        for f=1:Dclone.ntrials
            Dclone(spineind,cind(1,:),f)=real(squeeze(trialcspect(f,:,:)));
            Dclone(spineind,cind(2,:),f)=imag(squeeze(trialcspect(f,:,:)));
        end
    else
        Dclone=clone(DAll,cname,[size(DAll,1) Nclonesamples 1]);
        Dclone(spineind,cind(1,:),1)=real(cspect.crsspctrm);
        Dclone(spineind,cind(2,:),1)=imag(cspect.crsspctrm);
    end
    %% load optical scan stl

    [F, V]=stlread(backshape);
    sub=[];
    sub.vertices = V;
    sub.faces = F;
    sub_red = reducepatch(sub,0.5);

    %put this in fieldtrip structure
    subject=[];
    subject.pos=sub_red.vertices;
    subject.tri=sub_red.faces;

    %% Now set up a source space  cylinder

    allF=[];
    allR2=allF;
    allVE=allF;
    allmaxFstat=allF;

    %arbvals=strvcat('singleshell','George','localspheres','singlesphere','infinite');

    cyl_source=ft_read_headshape(cyl);

    arbvals=strvcat('infinite');


    for arbind=1:size(arbvals,1) %loop through volume conductors

        res=10; %cylinder grid resolution

        src=make_spine_grid_MES_cylinder(sub_red,cyl_source,res);

        %% now set up torso model, this gives lead field
        % mapping each point (and orientation) in source space to sensors

        Lf=[];

        switch deblank(arbvals(arbind,:))

            case {'singleshell','singlesphere','localspheres','infinite'}

                cfg = [];
                cfg.method = deblank(arbvals(arbind,:));
                cfg.grad=grad;
                cfg.conductivity = 1;
                cfg.unit='mm';
                [bmodel] = ft_prepare_headmodel(cfg,subject);


                cfg = [];
                cfg.headmodel = bmodel;
                cfg.sourcemodel = src;
                cfg.reducerank= 'no';
                [sourcemodel] = ft_prepare_leadfield(cfg, grad);

                %unravel lf

                startidx=1;
                for k=1:size(src.pos,1) %for each source point
                    Lf(:,startidx:startidx+2)=sourcemodel.leadfield{k};
                    startidx=startidx+3;
                end

            case 'George'

                [Lf,grad,src]=bem_leadfields_neckcast(subject,sub_fids,grad,src);
                %Lf dimensions are nchans x 3nsourcepts
            otherwise
                error('no forward model defined')

        end % switch

        D=spm_opm_attach_leadfield(grad.label,Lf,1,Dclone);

        if length(src.inside)<length(src.pos)
            warning('Not all sources in the space')
        end

        %% now run inversion
        origbad=D.badchannels;
        D = badchannels(D,1:D.nchannels, 1); %% set all channels to bad
        D = badchannels(D,spineind, 0); %% set spine channels to good
        D = badchannels(D,origbad, 1); %% set original bad channels to bad

        Nchans=length(spineind);

        D.inv{1}.inverse=[];
        D.inv{1}.inverse.type=invtype;

        D.inv{1}.inverse.woi=[cind(1,1)/D.fsample cind(1,end)/D.fsample].*1000; % in ms
        D.inv{1}.inverse.lpf=0;
        D.inv{1}.inverse.hpf=0;
        D.inv{1}.inverse.Han=0;

        if FIXMODES
            D.inv{1}.inverse.Nm=Nchans;
            D.inv{1}.inverse.A=eye(Nchans); %% force use of all channels
        else
            D.inv{1}.inverse.Nm=[];
            D.inv{1}.inverse.A=[]; 
        end

        D.inv{1}.inverse.no_temporal_filter=1;
        D.inv{1}.inverse.complexind=cind;
        D.inv{1}.forward.modality='MEG';
        D.inv{1}.forward.sensors=grad;
        D.inv{1}.forward.siunits=1; %check implications of this--> units
        D.val=1;
        D.inv{1}.inverse.Nt=[] ;


        Dout=TEST_spm_eeg_invert_classic_volumetric(D,1);

        %% now we have a mapping between source and sensor space given by M

        Ic=Dout.inv{1}.inverse.Ic{1}; %% channels
        F=Dout.inv{1}.inverse.F; %% free energy
        R2=Dout.inv{1}.inverse.R2; %% variance explained
        VE=Dout.inv{1}.inverse.VE; %% variance projected
        U=Dout.inv{1}.inverse.U{1}; %% reduction of data
        M=real(Dout.inv{1}.inverse.M*U); %% mapping from sensors to sources


        allF(arbind)=F;
        allR2(arbind)=R2;
        allVE(arbind)=VE;
    end % for arbind

    %% Now get per trial current density estimates in source space

    %frequencies of interest
    minf=freqroi(1);
    maxf=freqroi(2);
    fcind=intersect(find(fdat.freq>=minf),find(fdat.freq<=maxf));
    usefreq=fdat.freq(fcind);

    bdata=squeeze(fdat_brain.fourierspctrm(:,bind,fcind));
    spdata=squeeze(fdat.fourierspctrm(:,spind,fcind));
    emdata=squeeze(fdat_brain.fourierspctrm(:,fbemgind,fcind));
    ntapers=size(fdat.fourierspctrm,1);

    Nx=length(src.pos);
    J=zeros(Nx*3,length(fcind),ntapers);

    Jcov=zeros(Nx*3,Nx*3);
    for f=1:ntapers
        Jtrial=M*squeeze(spdata(f,:,fcind));
        Jcov=Jcov+cov(Jtrial');
        J(:,:,f)=Jtrial;
    end
    Jcov=Jcov./ntapers;

    xind=1:3:Nx*3; %indices for x y and z oriented sources
    yind=2:3:Nx*3;
    zind=3:3:Nx*3;

    %% sum power over all axes
    magJ=sqrt((diag(Jcov(xind,xind)+diag(Jcov(yind,yind)+diag(Jcov(zind,zind))))));
    %[maxJ,peakind]=max(magJ);

    %% get optimal orientation based on whole cord
    Jv=sqrt([sum(diag(Jcov(xind,xind))) sum(diag(Jcov(yind,yind))) sum(diag(Jcov(zind,zind)))]);
    Jv=Jv./sqrt(dot(Jv,Jv)); %this is optimal orientation

    if ~isempty(FIXORIENT)
        Jv=FIXORIENT;
        warning('Fixing orientation');
    end

    Jorient=J(xind,:,:).*Jv(1)+J(yind,:,:).*Jv(2)+J(zind,:,:).*Jv(3); %source data for opt orientation or fixorient
    Jocov=zeros(length(xind),length(xind));

    for f=1:ntapers
        Jtrial=Jorient(:,:,f);
        Jocov=Jocov+cov(Jtrial');
    end

    Jocov=Jocov./ntapers;
    [maxJ,peakindv]=max(diag(Jocov));

    magx=sqrt(diag(Jcov(xind,xind)));
    magy=sqrt(diag(Jcov(yind,yind)));
    magz=sqrt(diag(Jcov(zind,zind)));

    magopt=sqrt(diag(Jocov));

    J=J-mean(J,3);
    %Jorient=Jorient-mean(Jorient,3); %% along one axis

    %% Ybrain should be data from optimal linear mixture of channels to get emg coherence
    % Yibrain should be data from channel mixture orthogonal to this

    Ybrain=[];Yibrain=[];
    for f=1:size(fdat_brain.fourierspctrm,1)
        Ybrain(f,:)=Ur(:,1)'*squeeze(fdat_brain.fourierspctrm(f,bind,:));
        Yibrain(f,:)=Uinv(:,1)'*squeeze(fdat_brain.fourierspctrm(f,bind,:));
    end

    %% X is always cord i.e. Jorient (or J which is 3*longer)

    if COMB3AXES
        maxk=Nx;
        useJ=Jorient; %% optimal axis across whole cord
    else
        maxk=Nx*Nz*3;  %% 3 orientations evaluated separately
        useJ=J;
    end

    X=useJ;
    %% Y can either be emg or brain
    Y=[];
    if findstr(whatstr,'brain')
        Y=[real(Ybrain) imag(Ybrain)];
    end

    if findstr(whatstr,'emg')
        Y=[real(emdata) imag(emdata)];
    end
    %% keep copies of brain and emg signal for use later

    Ybrain2=[real(Ybrain) imag(Ybrain)];
    Zemg=[real(emdata) imag(emdata)];
    if isempty(Y)
        error('Y not defined');
    end

    %% linear regression: Y and ortho brain signal

    Fstatlin=zeros(length(usefreq),maxk);
    Yconfound=[real(Yibrain) imag(Yibrain)]; %% comes from orthogonal brain channel mixture

    % removing orthogonal brain signal
    Yr1=zeros(size(Y));
    rpcind=[cind cind];
    Forthri=zeros(size(cind,2),1);
    for f=1:size(cind,2)*2
        [B,BINT,Yr1(:,f),RINT,STATS] = regress(Y(:,f),[Yconfound(:,rpcind(1,f)) Yconfound(:,rpcind(2,f)) ones(size(Yconfound,1),1)]);
        Forthri(f)=STATS(2);
    end

    Forth_real=Forthri(cind(1,:)); %relationship between confound and signal
    Forth_imag=Forthri(cind(2,:)); %% separate real and imag

    if RMORTHBRAIN %Y becomes residuals if remove ortho brain (not for emg - cord)
        stvar=var(Y);
        Y=Yr1;
        fprintf('\n Mean %3.2f percent variance removed by orthogonal channels\n',100*(1-mean(var(Y)./stvar)))
    end

    %% get linear emg-brain
    rYbrain=zeros(size(Ybrain2));
    Fstatlin_emgbrainri=zeros(size(cind,2),1);
    pvals_emgbrainri=zeros(size(cind,2),1);


    for f=1:size(cind,2)*2

        %% Ybrain2 is copy of original brain signal, Zemg is copy of original EMG. Both are complex.
        [B,BINT,R,RINT,STATS] = regress(Ybrain2(:,f),[Zemg(:,rpcind(1,f)) Zemg(:,rpcind(2,f)) ones(size(Zemg,1),1)]);
        rYbrain(:,f)=R;
        Fstatlin_emgbrainri(f)=STATS(2);
        pvals_emgbrainri(f)=STATS(3);
        % testing at each freq so will have to deal with MC

    end % for f

    Fstatlin_emgbrain_real=Fstatlin_emgbrainri(cind(1,:));
    Fstatlin_emgbrain_imag=Fstatlin_emgbrainri(cind(2,:)); %% separate real and imag

    pvals_emgbrain_real=pvals_emgbrainri(cind(1,:));
    pvals_emgbrain_imag=pvals_emgbrainri(cind(2,:));

    % work out the threshold using FDR
    [~, crit_p_real, ~, adj_p_real]=fdr_bh(pvals_emgbrain_real,0.01,'dep','yes');
    [~, crit_p_imag, ~, adj_p_imag]=fdr_bh(pvals_emgbrain_imag,0.01,'dep','yes');
    realSig=find(adj_p_real<0.05);
    imagSig=find(adj_p_imag<0.05);

    figure;
    % Plot F statistics for real and imaginary parts of brain-EMG
    subplot(2, 1, 1);
    plot(usefreq, Fstatlin_emgbrain_real, 'color',col3,'LineWidth',3,'LineStyle',':'); hold on
    plot(usefreq, Fstatlin_emgbrain_imag, 'color',col4, 'LineWidth', 3);
    plot(usefreq(realSig), ones(1,length(realSig))*10,'*','MarkerSize',8,'color',col3)
    plot(usefreq(imagSig), ones(1,length(imagSig))*12,'*','MarkerSize',8,'color',col4)
    ylabel('F Statistic');
    xlabel('Frequency (Hz)');
    legend('Brain-EMG Real', 'Brain-EMG Imaginary');
    ax = gca;
    ax.FontSize = 14;
    box off
    xlim([0 40])
    % Plot Coherence between brain and EMG
    subplot(2, 1, 2);
    plot(brainmix_emg_coh.freq, brainmix_emg_coh.cohspctrm, 'LineWidth', 3,'color',col5);
    ylabel('Coherence');
    ax = gca;
    ax.FontSize = 14;
    box off
    xlim([0 40])
    sgtitle('Linear Brain-EMG', 'FontSize', 16, 'FontWeight', 'bold');
    set(gcf, 'Position', [100, 100, 600, 800]);  % Adjust figure size and position

    linvaremgbrain_rmvd=1-var(rYbrain(:))/var(Ybrain2(:)); %residuals as prop of orig brain signal

    %% now non-linear brain-emg cva (using residuals from linear regression)
    fprintf('\n Taking abs for brain residuals and emg') % THIS BIT OF CODE HARDCODED TO ABS
    rYbrain_abs=abs([rYbrain(:,cind(1,:))+i*rYbrain(:,cind(2,:))]); %% put it back into complex numbers then abs
    Zemg_abs=abs([Zemg(:,cind(1,:))+i*Zemg(:,cind(2,:))]);

    %mean centre
    rYbrain_abs=rYbrain_abs-mean(rYbrain_abs);
    Zemg_abs=Zemg_abs-mean(Zemg_abs);

    %CVA for abs brain-emg
    CVAemgbrain_abs=spm_cva(rYbrain_abs,Zemg_abs);
    fprintf('EMG brain CVA chi2=%.2f p=%.3f\n',CVAemgbrain_abs.chi(1),CVAemgbrain_abs.p(1))


    if SHUFFLE
        warning('SHUFFLING DATA !')
        rYbrain_abs=rYbrain_abs(randperm(size(rYbrain,1)),:);
        rYbrain_abs=rYbrain_abs-mean(rYbrain_abs);
        CVA_shuf_emgbrain=spm_cva(rYbrain_abs,Zemg_abs);
        fprintf('CVA shuf data emg brain chi2=%.2f p=%.3f \n', CVA_shuf_emgbrain.chi(1),CVA_shuf_emgbrain.p(1))
    end

    Ncan_emgbrain=max(find(CVAemgbrain_abs.p<0.05));
    normV_emgbrain=(cov(CVAemgbrain_abs.Y))*CVAemgbrain_abs.V(:,1:Ncan_emgbrain)*inv(cov(CVAemgbrain_abs.v(:,1:Ncan_emgbrain)));
    normW_emgbrain=(cov(CVAemgbrain_abs.X))*CVAemgbrain_abs.W(:,1:Ncan_emgbrain)*inv(cov(CVAemgbrain_abs.w(:,1:Ncan_emgbrain)));

    figure
    plot(usefreq,abs(normV_emgbrain(:,1)),'LineWidth',3,'color',col1);
    xlabel('Frequency (Hz)','FontSize',20)
    ylabel('Brain')
    box off
    ax = gca;
    ax.FontSize = 20;
    hold on
    yyaxis right
    plot(usefreq,abs(normW_emgbrain(:,1)),'LineWidth',3,'color', col2);
    ylabel('EMG')
    legend({'Brain', 'EMG'})

    %% now cord and Y linear regression, scannning up cord

    Y=Y-mean(Y);
    Fstatlinri=zeros(size(cind,2),maxk);
    pvalsri=zeros(size(cind,2),maxk);
    rYcord=zeros([size(Y),maxk]); %% need specific residual at each point on cord
    for k=1:maxk
        X=[real(squeeze(useJ(k,:,:))') imag(squeeze(useJ(k,:,:))')]; % cord
        X=X-mean(X);

        for f=1:size(cind,2)*2 %% want to deal with real and imag part of Y

            [B,BINT,R,RINT,STATS] = regress(Y(:,f),[X(:,rpcind(1,f)) X(:,rpcind(2,f)) ones(size(X,1),1)]);
            rYcord(:,f,k)=R; %% THIS WAS WRONG AS rY changed at each location on cord...

            % note also that rYcord has both real and imag components that had linear relationshop removed

            Fstatlinri(f,k)=STATS(2);
            rsqlin(f,k)=STATS(1);
            pvalsri(f,k)=STATS(3);
        end % for f

    end % for k

    Fstatlin_real=Fstatlinri(cind(1,:),:); %% separate real and imag
    Fstatlin_imag=Fstatlinri(cind(2,:),:);
    pvals_real=pvalsri(cind(1,:),:);
    pvals_imag=pvalsri(cind(2,:),:);

 
    Y2=Y; %save a copy to use to calc variance later on; Y is redefined

    switch whatstr %% flags  RMORTHBRAIN and REMOVELIN now take care of options

        case 'brainopt'
            warning('not checked')
            Ycord=rYcord(:,cind(1,:),:) +i*rYcord(:,cind(2,:),:);
            X0=[]; %% already removed by RMORTH flag

        case 'brainopt+abs'
            %% Y will already have had linear projection from orth brain removed (if selected)
            % and also linear
            Ycord=abs(rYcord(:,cind(1,:),:) +i*rYcord(:,cind(2,:),:));

            X0=abs(Yibrain); % % abs of orth brain : use confound here as non-linear portion may influence

        case 'orthbrain+brainopt'
            warning('not checked')
            if RMORTHBRAIN
                error('orthbrain already removed as confound');
            end
            X0=Y; %% null space is optimal brain signal
            Y=[real(Yibrain) imag(Yibrain)]; %% signal space is orthogonal (to optimal) brain signal

        case 'emg'
            warning('not checked')
            Ycord=rYcord(:,cind(1,:),:) +i*rYcord(:,cind(2,:),:);
            X0=[];
        case 'emg+abs'
            warning('not checked')
            %Yc=abs([Y(:,cind(1,:))+i*Y(:,cind(2,:))]); %% put it back into complex numbers then abs
            %Y=abs(Yc);
            Ycord=abs(rYcord(:,cind(1,:),:) +i*rYcord(:,cind(2,:),:));
            X0=abs(Yibrain); % abs of orth brain : use confound here as non-linear portion may influence

        otherwise
            error('not defined')
    end % switch

    %% CVA with linear part removed
    X0=X0-mean(X0);
    chi2=0;
    chi2r=0;

    pval=[];
    allCVA=[];

    for k=1:maxk

        X=squeeze(useJ(k,:,:))'; % cord data

        if contains(whatstr,'abs')
            X=abs(X);
        else
            X=[real(X) imag(X)];
        end

        %% here is where we need to account for/ quantify linear effect from part of cord (i.e. just before non lin stage)
        Y=squeeze(Ycord(:,:,k)); %% this will have had linear portion removed if required

        Y=Y-mean(Y);
        X=X-mean(X);

        CVA=spm_cva(Y,X,X0);
        chi2(k)=CVA.chi(1);
        pval(k)=CVA.p(1);
        cvar2(k)=CVA.r(1).^2;
        allCVA{k}=CVA;
    end % for k


    if COMB3AXES

       [~,peakind]=max(magopt); %% peak power
       clear maxk

       if contains(filenames(cnd,:),'124') %need second max as well for 123
        [~,maxidx2]=maxk(magopt,2);
       end

       figure;
       plot_func_spine_dat(subject,src,magopt,grad,col1)
       title('power')
       hold on
       plot3(src.pos(peakind,1), src.pos(peakind,2),src.pos(peakind,3),'ro','MarkerSize',10,'MarkerFaceColor','r')
      
       if contains(filenames(cnd,:),'124')
           plot3(src.pos(maxidx2(2),1), src.pos(maxidx2(2),2),src.pos(maxidx2(2),3),'ro','MarkerSize',10,'MarkerFaceColor','r')
           peakind=maxidx2(2);
       end

        %calc variance removed for point of max power
        rYMax=squeeze(rYcord(:,:,peakind));
        linvar_rmvd= 1-var(rYMax(:))/var(Y2(:)); %save this for point of max pwr. Y2 is a copy of Y from lin reg (is overwritten)

        CVA=allCVA{peakind}; %cva for sourcepoint where power is max
        fprintf('CVA at point of peak power chi2=%.2f p=%.3f \n', CVA.chi(1),CVA.p(1))
        fprintf('Source recon VE is %.2f R2 is %.2f\n', allVE,allR2)

        Ncan=max(find(CVA.p<0.05));

        normV=cov(CVA.Y)*CVA.V(:,1:Ncan)*inv(cov(CVA.v(:,1:Ncan))); %
        normW=cov(CVA.X)*CVA.W(:,1:Ncan)*inv(cov(CVA.w(:,1:Ncan))); %

        %shuffle just for point of max pwr

        if SHUFFLE
            warning('SHUFFLING DATA !')
            Y=squeeze(Ycord(:,:,peakind)); %% this will have had linear portion removed if required
            Y=Y(randperm(size(Y,1)),:);
            X=squeeze(useJ(peakind,:,:))'; % cord data
            Y=Y-mean(Y);
            X=X-mean(X);
            CVA_shuf=spm_cva(Y,X,X0);
            fprintf('CVA shuf data at point of peak power chi2=%.2f p=%.3f \n', CVA_shuf.chi(1),CVA_shuf.p(1))
        end

        cvaname=[savepath filesep sprintf('%s_%s.mat',pstr,whatstr)];
        save(cvaname,'allCVA','peakind','magopt','usefreq','CVAemgbrain_abs','Fstatlinri','Fstatlin_emgbrainri','linvar_rmvd','linvaremgbrain_rmvd');

        %% plot nonlinear

        cols=colormap(brewermap([],"Dark2"));
        col1=cols(1,:); col2=cols(6,:);
        figure;
        plot(usefreq,abs(normV(:,1)),'LineWidth',3,'color',col1);
        if contains(whatstr,'emg')
            leg1='EMG';
            ylabel('EMG')
        else
            leg1='Brain';
            ylabel('Brain')
        end
        xlabel('Frequency (Hz)','FontSize',20)
        ax = gca;
        ax.FontSize = 20;
        yyaxis right
        ylabel('Spinal cord')
        plot(usefreq,abs(normW(:,1)),'LineWidth',3,'color', col2);
        legend({leg1,'Spinal Cord'})
        box off
        
        figure
        subplot(2,1,1);
        h=plot(usefreq,Fstatlin_real(1:length(usefreq),peakind),usefreq,Forth_real);
        set(h(1),'Linewidth',4);
        title('real')

        subplot(2,1,2);
        h=plot(usefreq,Fstatlin_imag(1:length(usefreq),peakind),usefreq,Forth_imag);
        set(h(1),'Linewidth',4);
        title('imag')
        if RMORTHBRAIN
            legend('Linear post removal','Brain orth removed');
        else
            legend('Linear','Brain orth not removed');
        end
        %savename=['D:\figsfortalk\nonlin_left_spineemg.png'];
        %exportgraphics(gcf, savename, 'Resolution', 600)

        %% coherence for visualization
        SC_freqdat=squeeze(useJ(peakind,:,:))'; % cord, complex

        Y=[];
        if findstr(whatstr,'brain')
            Y=Ybrain;
        end

        if findstr(whatstr,'emg')
            Y=emdata;
        end

        [Coh_ft] = cohxy(SC_freqdat,Y,fdat);
        
        p_real=pvals_real(:,peakind);
        p_imag=pvals_imag(:,peakind);

        % work out the threshold using FDR
        [~, crit_p_r, ~, adj_p_real]=fdr_bh(p_real,0.01,'dep','yes');
        [~, crit_p_i, ~, adj_p_imag]=fdr_bh(p_imag,0.01,'dep','yes');
        realSig=find(adj_p_real<0.05);
        imagSig=find(adj_p_imag<0.05);

        figure; %plot coherence and f stats in subplot
        subplot(2, 1, 1);
        plot(usefreq, Fstatlin_real(1:length(usefreq), peakind), 'LineWidth', 3,'color',col3); hold on
        plot(usefreq, Fstatlin_imag(1:length(usefreq), peakind), 'LineWidth', 3,'color',col4,'LineStyle',':');
        plot(usefreq(realSig), ones(1,length(realSig))*400,'*','MarkerSize',8,'color',col3)
        plot(usefreq(imagSig), ones(1,length(imagSig))*450,'*','MarkerSize',8,'color',col4)
        legend({'Real', 'Imaginary'})
        xlabel('Frequency (Hz)');
        ylabel('F statistic');
        box off
        ax = gca;
        ax.FontSize = 12;

        subplot(2, 1, 2);
        plot(usefreq, squeeze(Coh_ft.cohspctrm(1, 2, :)), 'LineWidth', 3,'color',col5);
        xlabel('Frequency (Hz)');
        ylabel('Coherence');
        set(gcf, 'Position', [100, 100, 600, 800]); % Set figure size and position
        box off
        ax = gca;
        ax.FontSize = 12;
        sgtitle('Linear Cord-Y', 'FontSize', 16, 'FontWeight', 'bold');
        % savename=['D:\figsfortalk\linear_right_brainspine.png'];


    end % if COMB3AXES

end % for filenames



