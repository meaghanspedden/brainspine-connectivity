%% SPM inversion - new neck cast - hand contraction July 2023

%for subject 123
clear all;
close all;

mydir='D:\MSST001\';

%filenames for concatenated runs for right hand and left hand
filenames=strvcat([mydir 'sub-OP00124\ses-001\meg\pprhandoe1000msfdfflo45hi5ds_sub-OP00124_ses-001_task-staticright_run-001_meg'],...
        [mydir '\sub-OP00124\ses-001\meg\pplhandoe1000msfdfflo45hi5ds_sub-OP00124_ses-001_task-staticleft_run-001_meg']);

posfile=[mydir 'sub-OP00124\ses-001\meg\ds_sub-OP00124_ses-001_positions.tsv'];
backshape='D:\OP00123_experiment\cast space\OP00123.stl'; %stl in same space as sensors
cyl='D:\OP00123_experiment\cast space\cylinder_good_space.stl'; %cylinder for source space

dpath=[mydir 'sub-OP00124\ses-001\meg']; %data path
savepath=[mydir '\Coh_results00123']; %save path

if ~exist(savepath,'dir')
    mkdir(savepath)
end

cd(dpath)

brainchan_labels={'19','DG', 'OH', 'A1','1B', 'A9','JS','A6','DJ','MY','DS','OK','MI','17'};
badchans={'ML-X','ML-Y','ML-Z','K4-Z','K4-X'}; %%

%%
addpath D:\torso_tools
addpath D:\spm12

spm('defaults','EEG');
addpath(genpath('D:\brainspineconnectivity'))

%% analysis options
SHUFFLE=0;
allcanfilenames=[];

%whatstr='emg+abs';
whatstr='brainopt+abs';
%whatstr='orthbrain+brainopt';
%whatstr='emg';

FIXORIENT=[]; %
REMOVELIN=1;
RMORTHBRAIN=contains(whatstr,'brain'); %% don't remove brain if dealing with cord-muscle
FIXMODES=0;
CLONETRIALS=1; % use all trial data (1) rather than average (0)


freqroi=[5 35];

%invtype='EBBr95';
invtype='IID';
%invtype='EBB';

IMAGECROSS=0;
COMB3AXES=1;



sub_fids= [-201 238 123; %r shoulder
    -196 -207 157; %l shoulder
    -221 8.92 219.038  ];

for cnd=2 %:size(filenames,1),

    DAll=spm_eeg_load(filenames(cnd,:));
    [a1,b1,c1]=fileparts(DAll.fullfile);
    pstr=b1(1:10); %%get upper/lower l/r info from filename

    cname=[a1 filesep 'clone_b1' b1 c1 ];


    D=DAll; %

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
    cd      D:\fieldtrip-master\fieldtrip-master\private

    seedperm=0;
    cohbrain=coh_meaghan(D,D.chanlabels(brainind),D.chanlabels(emgind(1)),[min(freqroi) max(freqroi)]);
    [cspect_brain,fdat_brain,trialcspect_brain]=xspectrum_meaghan(D,D.chanlabels(brainind),D.chanlabels(emgind(1)),[min(freqroi) max(freqroi)],seedperm);
    [cspect,fdat,trialcspect]=xspectrum_meaghan(D,D.chanlabels(spineind),D.chanlabels(emgind(1)),[min(freqroi) max(freqroi)],seedperm);

    fbemgind=find(contains(fdat_brain.label,D.chanlabels(emgind(1))));
    fbrainind=setdiff(1:length(fdat_brain.label),fbemgind);

    %% plot some metrics of brain-emg interaction
    nchansBrain=length(cspect_brain.labelcmb); %it was getting confused because csd mat is square

    figure; subplot(4,1,1);
    plot(cspect_brain.freq,real(cspect_brain.crsspctrm(1:nchansBrain,:))');
    title('real cross')
    subplot(4,1,2);
    plot(cspect_brain.freq,imag(cspect_brain.crsspctrm(1:nchansBrain,:))');
    title('imag cross')
    subplot(4,1,3)
    plot(cspect_brain.freq,abs(cspect_brain.crsspctrm(1:nchansBrain,:))');
    title('abs cross')
    subplot(4,1,4)
    plot(cohbrain.freq,cohbrain.cohspctrm(1:nchansBrain,:)');
    title('coh')
    %% get dominant spatial component that explains the brain-emg cross spectrum
    r1=cspect_brain.crsspctrm;
    r1=r1-mean(r1);

    [Ur,S,V]=svd(real(r1)*real(r1)'); %% could be an issue if peak coherence is imag.
    varexp=cumsum(diag(S))./sum(diag(S));
    Nrcomp=min(find(varexp>0.99));
    %Ur=Ur(:,1:Nrcomp);
    [Ui,Si,Vi]=svd(imag(r1)*imag(r1)');

    Uinv=pinv(Ur); %% the spatial mixture of brain chans orthogonal to the ideal mixture for EMG-brain cross-spectrum
    CHECKMIX=1;
    %     if CHECKMIX
    %         figure;
    %         %% can base it on real or imag parts, very similar
    %         plot(1:length(S),Ur(:,1),1:length(S),Ui(:,1)) %% linear mixture of channels to give real or imag comp
    %
    %         replace.U=Ur(:,1); %% Ur or Uinv
    %         replace.ind=fbrainind;
    %         replace.label='G2-19-Y';
    %
    %         %% make a replace one of the chans with optimcal linear mixture and check it is larger/smaller
    %         [cspect_brainmix,fdat_brainmix,trialcspect_brainmix]=xspectrum_meaghan(D,D.chanlabels(brainind),D.chanlabels(emgind(1)),[min(freqroi) max(freqroi)],seedperm,[],replace);
    %         figure;
    %         plot(cspect_brain.freq,abs(cspect_brain.crsspctrm),'r');
    %         hold on;
    %         plot(cspect_brain.freq,abs(cspect_brainmix.crsspctrm),'go');
    %         ind=find(contains(cspect_brainmix.labelcmb,replace.label));
    %         plot(cspect_brain.freq,abs(cspect_brainmix.crsspctrm),'go');
    %         plot(cspect_brain.freq,abs(cspect_brainmix.crsspctrm(ind,:)),'k*');
    %         title('abs cross')
    %     end % checkmix

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
    %% plot back and sensors on back

    [F, V]=stlread(backshape);
    %     [a1,b1,c1]=fileparts(backshape);
    sub=[];
    sub.vertices = V;
    sub.faces = F;
    sub_red = reducepatch(sub,0.5);

    %put this in fieldtrip structure
    subject=[];
    subject.pos=sub_red.vertices;
    subject.tri=sub_red.faces;


    %% Now set up a source space either flat plane or cylinder

    allF=[];
    allR2=allF;
    allVE=allF;
    allmaxFstat=allF;


    %arbvals=strvcat('singleshell','George','localspheres','singlesphere','infinite');

    cyl_source=ft_read_headshape(cyl);

    arbvals=strvcat('infinite');


    for arbind=1:size(arbvals,1)


        count=0;

        res=10;

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
                [bmodel,cfg] = ft_prepare_headmodel(cfg,subject);


                cfg = [];
                cfg.headmodel = bmodel;
                cfg.sourcemodel = src;
                cfg.reducerank= 'no';
                [sourcemodel] = ft_prepare_leadfield(cfg, grad); %% takes a long time- probably due to dense back mesh- but concise output

                %put leadfields into mat similar to bem lf

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

        if ~isempty(badgrad) %have not tested georges
            warning('have not tested this works for BEM')

            badidx=find(contains(grad.label,badchans));
            Lf=Lf(setdiff(1:size(Lf,1),badidx),:); %% take out bad channels post-hoc
        end


        if length(src.inside)<length(src.pos)
            warning('Not all sources in the space')
        end

        %% GRB TEST
        %% now run inversion
        origbad=D.badchannels;
        D=Dclone;
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
            D.inv{1}.inverse.A=[]; % eye(Nchans); %% force use of all channels
        end

        D.inv{1}.inverse.no_temporal_filter=1;
        D.inv{1}.inverse.complexind=cind;
        D.inv{1}.forward.modality='MEG';
        D.inv{1}.forward.sensors=grad;
        D.inv{1}.forward.siunits=1; %I think everyting needs to be in
        %m for this to work

        D.val=1;
        D.inv{1}.inverse.Nt=[] ;


        Dout=spm_eeg_invert_classic_volumetric(D,1,Lf);

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


    % figure
    % bar(allF); xticklabels(arbvals)
    % ylabel('FE')
    % figure
    % bar(allR2); xticklabels(arbvals)
    % ylabel('R2')
    % figure
    % bar(allVE); xticklabels(arbvals)
    % ylabel('VE')

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

    %I think nx here is num sourcepts
    Nx=length(src.pos);
    J=zeros(Nx*3,length(fcind),ntapers);


    %magM=sqrt(diag((M*M')));
    %nM=M./repmat(magM,1,size(M,2)); % normalize M

    Jcov=zeros(Nx*3,Nx*3);
    for f=1:ntapers
        Jtrial=M*squeeze(spdata(f,:,fcind));
        Jcov=Jcov+cov(Jtrial');
        J(:,:,f)=Jtrial;
    end
    Jcov=Jcov./ntapers;


    xind=1:3:Nx*3;
    yind=2:3:Nx*3;
    zind=3:3:Nx*3;

    %% sum power over all axes
    magJ=sqrt((diag(Jcov(xind,xind)+diag(Jcov(yind,yind)+diag(Jcov(zind,zind))))));
    [maxJ,peakind]=max(magJ);

    %% get optimal orientation based on whole cord
    Jv=sqrt([sum(diag(Jcov(xind,xind))) sum(diag(Jcov(yind,yind))) sum(diag(Jcov(zind,zind)))]);
    Jv=Jv./sqrt(dot(Jv,Jv)); %this is optimal orient

    plotOptOri(Jv,src,subject)


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

    figure
    plot_func_spine_dat(subject,src,magopt,grad)
    title('Opt oriented sources')


    figure
    plot_func_spine_dat(subject,src,magx,grad)
    title('X oriented sources')
    %
    figure
    plot_func_spine_dat(subject,src,magz,grad)
    title('Z oriented sources')

    figure
    plot_func_spine_dat(subject,src,magy,grad)
    title('Y oriented sources')




    J=J-mean(J,3);
    %Jorient=Jorient-mean(Jorient,3); %% along one axis

    %% Ybrain should be data from optimal linear mixture of channels to get emg coherence
    %% Yibrain should be data from channel mixture orthogonal to this
    Ybrain=[];Yibrain=[];
    for f=1:size(fdat_brain.fourierspctrm,1),
        Ybrain(f,:)=Ur(:,1)'*squeeze(fdat_brain.fourierspctrm(f,bind,:));
        Yibrain(f,:)=Uinv(:,1)'*squeeze(fdat_brain.fourierspctrm(f,bind,:));
    end;



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
    end;


    if findstr(whatstr,'emg')
        Y=[real(emdata) imag(emdata)];
    end;
    %% keep copies of brain and emg signal for use later
    Ybrain2=[real(Ybrain) imag(Ybrain)]; %%
    Zemg=[real(emdata) imag(emdata)]; %% may need this later

    if isempty(Y),
        error('Y not defined');
    end;

    %% now optionally remove linear parts
    Fstatlin=zeros(length(usefreq),maxk);
    Yconfound=[real(Yibrain) imag(Yibrain)]; %% comes from orthogonal brain channel mixture

    %% removing orthogonal brain signal
    Yr1=zeros(size(Y));
    rpcind=[cind cind];
    Forthri=zeros(size(cind,2),1)
    for f=1:size(cind,2)*2,
        [B,BINT,Yr1(:,f),RINT,STATS] = regress(Y(:,f),[Yconfound(:,rpcind(1,f)) Yconfound(:,rpcind(2,f)) ones(size(Yconfound,1),1)]);
        Forthri(f)=STATS(2);
    end;
    Forth_real=Forthri(cind(1,:));
    Forth_imag=Forthri(cind(2,:)); %% separate real and imag

    if RMORTHBRAIN,
        stvar=var(Y);
        Y=Yr1;
        fprintf('\n Mean %3.2f percent variance removed by orthogonal channels\n',100*(1-mean(var(Y)./stvar)))
    end;
    %% get linear emg-brain
    rYbrain=zeros(size(Ybrain2));

    Fstatlin_emgbrainri=zeros(size(cind,2),1);
    for f=1:size(cind,2)*2,
        %% Ybrain2 is copy of original brain signal, Zemg is copy of original EMG. Both are complex.
        [B,BINT,R,RINT,STATS] = regress(Ybrain2(:,f),[Zemg(:,rpcind(1,f)) Zemg(:,rpcind(2,f)) ones(size(Zemg,1),1)]);
        rYbrain(:,f)=R;

        Fstatlin_emgbrainri(f)=STATS(2);

    end; % for f

    Fstatlin_emgbrain_real=Fstatlin_emgbrainri(cind(1,:))
    Fstatlin_emgbrain_imag=Fstatlin_emgbrainri(cind(2,:)); %% separate real and imag

    figure;
    plot(usefreq,Fstatlin_emgbrain_real,usefreq,Fstatlin_emgbrain_imag);
    legend('Brain-emg real','Brain-emg imag');
    title('linear');
    linvaremgbrain_rmvd=1-var(rYbrain(:))/var(Ybrain2(:));

    % now non-linear brain emg
    fprintf('\n Taking abs for brain residuals and emg') % THIS BIT OF CODE HARDCODED TO ABS
    rYbrain_abs=abs([rYbrain(:,cind(1,:))+i*rYbrain(:,cind(2,:))]); %% put it back into complex numbers then abs
    Zemg_abs=abs([Zemg(:,cind(1,:))+i*Zemg(:,cind(2,:))]);
    rYbrain_abs=rYbrain_abs-mean(rYbrain_abs);
    Zemg_abs=Zemg_abs-mean(Zemg_abs);
    CVAemgbrain_abs=spm_cva(rYbrain_abs,Zemg_abs);
    
     Ncan_emgbrain=max(find(CVAemgbrain_abs.p<0.05));
     normV_emgbrain=(cov(CVAemgbrain_abs.Y))*CVAemgbrain_abs.V(:,1:Ncan_emgbrain);
     normW_emgbrain=(cov(CVAemgbrain_abs.X))*CVAemgbrain_abs.W(:,1:Ncan_emgbrain);

        cols=colormap(brewermap([],"Dark2"));
        col1=cols(1,:); col2=cols(6,:);
        subplot(211)
        plot(usefreq,abs(normV_emgbrain(:,1)),'LineWidth',3,'color',col1);
        title('Brain','FontSize',20)
    
        xlabel('Frequency (Hz)','FontSize',20)
        box off
        ax = gca;
        ax.FontSize = 20;

        subplot(212)
        plot(usefreq,abs(normW_emgbrain(:,1)),'LineWidth',3,'color', col2);
        title('EMG','FontSize',20)
        box off
        ax = gca;
        ax.FontSize = 20;
        xlabel('Frequency (Hz)','FontSize',20)

    %% now scanning up cord


    Y=Y-mean(Y);
    Fstatlinri=zeros(size(cind,2),maxk);
    rYcord=zeros([size(Y),maxk]); %% need specific residual at each point on cord
    for k=1:maxk,
        k
        X=[real(squeeze(useJ(k,:,:))') imag(squeeze(useJ(k,:,:))')]; % cord
        X=X-mean(X);

        for f=1:size(cind,2)*2, %% want to deal with real and imag part of Y

            [B,BINT,R,RINT,STATS] = regress(Y(:,f),[X(:,rpcind(1,f)) X(:,rpcind(2,f)) ones(size(X,1),1)]);
            rYcord(:,f,k)=R; %% THIS WAS WRONG AS rY changed at each location on cord...
            %% note also that rYcord has both real and imag components that had linear relationshop removed


            Fstatlinri(f,k)=STATS(2);
            rsqlin(f,k)=STATS(1);
        end; % for f

    end; % for k

    Fstatlin_real=Fstatlinri(cind(1,:),:); %% separate real and imag
    Fstatlin_imag=Fstatlinri(cind(2,:),:); %% separate real and imag
    figure;
    plot(usefreq,Fstatlin_real,usefreq,Fstatlin_imag,':');
    legend('linear real','linear imag');
    title('linear');


    switch whatstr %% flags  RMORTHBRAIN and REMOVELIN now take care of options

        case 'brainopt',
            warning('not checked')
            Ycord=rYcord(:,cind(1,:),:) +i*rYcord(:,cind(2,:),:);
            X0=[]; %% already removed by RMORTH flag
        case 'brainopt+abs',
            %% Y will already have had linear projection from orth brain removed (if selected)
            %% and also linear
            %Y=abs([Y(:,cind(1,:))+i*Y(:,cind(2,:))]); %% put it back into complex numbers then abs
            Ycord=abs(rYcord(:,cind(1,:),:) +i*rYcord(:,cind(2,:),:));

            X0=abs(Yibrain); % % abs of orth brain : use confound here as non-linear portion may influence

        case 'orthbrain+brainopt',
            warning('not checked')
            if RMORTHBRAIN,
                error('orthbrain already removed as confound');
            end;
            X0=Y; %% null space is optimal brain signal
            Y=[real(Yibrain) imag(Yibrain)]; %% signal space is orthogonal (to optimal) brain signal

        case 'emg',
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
    end; % switch


    %
    % if size(X)~=size(Y),
    %     error('X and Y don not match')
    % end;
    %
    % if ~isempty(X0)
    %     if size(X0)~=size(X),
    %         error('X and X0 do not match')
    %     end;
    % end;

    % Y=Y-mean(Y);
    % if SHUFFLE,
    %     warning('SHUFFLING DATA !')
    %     Y=Y(randperm(size(Y,1)),:);
    % end;

    X0=X0-mean(X0);
    chi2=0;
    chi2r=0;

    pval=[];
    allCVA=[];

    for k=1:maxk,
        k

        X=squeeze(useJ(k,:,:))'; % cord data


        if contains(whatstr,'abs'),
            X=abs(X);
        else
            X=[real(X) imag(X)];
        end;

        %% here is where we need to account for/ quantify linear effect from part of cord (i.e. just before non lin stage)
        Y=squeeze(Ycord(:,:,k)); %% this will have had linear portion removed if required
        Y=Y-mean(Y);
        X=X-mean(X);

        CVA=spm_cva(Y,X,X0);
        chi2(k)=CVA.chi(1);
        pval(k)=CVA.p(1);
        cvar2(k)=CVA.r(1).^2;
        allCVA{k}=CVA;
    end; % for k


    if COMB3AXES

        PWRMAX=1;

        clear maxk
        [maxval maxidx] = maxk(magopt,2);
        if PWRMAX
            [dum,peakind]=max(magopt); %% peak power
        else
            [dum,peakind]=max(chi2); %% peak stat
        end

    

        figure;
        plot_func_spine_dat(subject,src,magopt,grad)
        title('power')
        hold on
        plot3(src.pos(peakind,1), src.pos(peakind,2),src.pos(peakind,3),'ro','MarkerSize',12)
        plot3(src.pos(maxidx(2),1), src.pos(maxidx(2),2),src.pos(maxidx(2),3),'ro','MarkerSize',12)

        peakind=maxidx(2);

%         [dum,peakindchi]=max(chi2); %% peak stat
%         figure;
%         plot_func_spine_dat(subject,src,chi2,grad)
%         title('chi2')
%         hold on
%         plot3(src.pos(peakindchi,1), src.pos(peakindchi,2),src.pos(peakindchi,3),'ro','MarkerSize',12)
%         title('CVA chi2')

        CVA=allCVA{peakind}; %cva for sourcepoint where power is max

        Ncan=max(find(CVA.p<0.05));
        normV=(cov(CVA.Y))*CVA.V(:,1:Ncan);
        normW=(cov(CVA.X))*CVA.W(:,1:Ncan);
        cvaname=[savepath filesep sprintf('%s_%s.mat',pstr,whatstr)];
        %save(cvaname,'allCVA','peakind','magopt','usefreq','CVAemgbrain','Fstatlin','Fstatlin_emgbrain','linvar_rmvd','linvaremgbrain_rmvd');



        %% plot nonlinear

        cols=colormap(brewermap([],"Dark2"));
        col1=cols(1,:); col2=cols(6,:);
        subplot(211)
        plot(usefreq,abs(normV(:,1)),'LineWidth',3,'color',col1);
        if contains(whatstr,'emg')
            title('EMG','FontSize',20)
        else
            title('Brain','FontSize',20)
        end
        xlabel('Frequency (Hz)','FontSize',20)
        box off
        ax = gca;
        ax.FontSize = 20;

        subplot(212)
        plot(usefreq,abs(normW(:,1)),'LineWidth',3,'color', col2);
        title('Cord','FontSize',20)
        box off
        ax = gca;
        ax.FontSize = 20;
        xlabel('Frequency (Hz)','FontSize',20)


        %savename=['D:\figsfortalk\nonlin_left_spineemg.png'];

        %  exportgraphics(gcf, savename, 'Resolution', 600)

        figure;
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
        end;


        ylabel('Fstat')
        xlabel('freq');



        figure;
        plot(usefreq,Fstatlin_real(1:length(usefreq),peakind),'LineWidth',3,'color','k');
        hold on;
        plot(usefreq,Fstatlin_imag(1:length(usefreq),peakind),'LineWidth',3,'color','r');
        legend('real','imag')
        ylabel('Fstat','Fontsize',20)
        xlabel('Frequency (Hz)','FontSize',20)
        ax = gca;
        ax.FontSize = 20;
        %legend('Left contraction','Right contraction','FontSize',20)
        box off
        %ylim([0 2.5])

        savename=['D:\figsfortalk\linear_right_brainspine.png'];

        
    end; % if COMB3AXES

end; % for filenames


