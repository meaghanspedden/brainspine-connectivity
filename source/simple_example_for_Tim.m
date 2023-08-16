%% SPM inversion - new neck cast - hand contraction

clear all;
close all;
clc

mydir='D:\MSST001';
subjectID ='123'; %122 or %123
%whatstr='emg+abs';
whatstr='brainopt+abs';

%% paths
fieldtrippath='D:\fieldtrip-master\fieldtrip-master';
addpath D:\spm12
spm('defaults','EEG');
addpath(genpath('D:\brainspineconnectivity'))

%% get meta data

dataOut=getMetaData(subjectID, mydir);

filenames=dataOut.filenames;
posfile=dataOut.posfile;
backshape=dataOut.backshape;
cyl=dataOut.cyl;
brainchan_labels=dataOut.brainchan_labels;
badchans=dataOut.badchans;
dpath=dataOut.dpath;
savepath=dataOut.savepath;


if ~exist(savepath,'dir')
    mkdir(savepath)
end

cd(dpath)

%% analysis options

FIXORIENT=[]; %empty means it calculates optimal orientation
FIXMODES=0;
CLONETRIALS=1; % use all trial data (1) rather than average (0)

freqroi=[5 35];
invtype='IID';

cnd=1;

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
    cd([fieldtrippath,'\private']) %need to fix this I think it is fieldtrip and spm interacting
    seedperm=0;
    cohbrain=coh_meaghan(D,D.chanlabels(brainind),D.chanlabels(emgind(1)),[min(freqroi) max(freqroi)]);
    [cspect_brain,fdat_brain,trialcspect_brain]=xspectrum_meaghan(D,D.chanlabels(brainind),D.chanlabels(emgind(1)),[min(freqroi) max(freqroi)],seedperm);
    [cspect,fdat,trialcspect]=xspectrum_meaghan(D,D.chanlabels(spineind),D.chanlabels(emgind(1)),[min(freqroi) max(freqroi)],seedperm);

    fbemgind=find(contains(fdat_brain.label,D.chanlabels(emgind(1))));
    fbrainind=setdiff(1:length(fdat_brain.label),fbemgind);


    %% get dominant spatial component that explains the brain-emg cross spectrum

    r1=cspect_brain.crsspctrm;
    r1=r1-mean(r1);

    [Ur,S,~]=svd(real(r1)*real(r1)'); %% use real part
    [Ui,Si,Vi]=svd(imag(r1)*imag(r1)');
    varexp=cumsum(diag(S))./sum(diag(S));
    Nrcomp=min(find(varexp>0.99));
    Uinv=pinv(Ur); %% the spatial mixture of brain chans orthogonal to the ideal mixture for EMG-brain cross-spectrum


    %calculate coherence with ideal mixture-EMG
    [~, fdat_brainmix,~]=xspectrum_meaghan(D,D.chanlabels(brainind),D.chanlabels(emgind(1)),[min(freqroi) max(freqroi)],seedperm,[],[]);



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

    %% Now set up a source space cylinder

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

        %% set up torso model, this gives lead field
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
                error('This is not set up yet')
                [Lf,grad,src]=bem_leadfields_neckcast(subject,sub_fids,grad,src);
                %Lf dimensions are nchans x 3nsourcepts
            otherwise
                error('No forward model defined')

        end % switch

        D=TEST_spm_opm_attach_leadfield(grad.label,Lf,1,Dclone);

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


    %% get optimal orientation based on whole cord
    Jv=sqrt([sum(diag(Jcov(xind,xind))) sum(diag(Jcov(yind,yind))) sum(diag(Jcov(zind,zind)))]);
    Jv=Jv./sqrt(dot(Jv,Jv)); %this is optimal orientation

    if ~isempty(FIXORIENT) %if user specifies orientation use this instead
        Jv=FIXORIENT;
    end

    Jorient=J(xind,:,:).*Jv(1)+J(yind,:,:).*Jv(2)+J(zind,:,:).*Jv(3); %source data for opt orientation or fixorient
    Jocov=zeros(length(xind),length(xind));

    for f=1:ntapers
        Jtrial=Jorient(:,:,f);
        Jocov=Jocov+cov(Jtrial');
    end

    Jocov=Jocov./ntapers;

    %     magx=sqrt(diag(Jcov(xind,xind)));
    %     magy=sqrt(diag(Jcov(yind,yind)));
    %     magz=sqrt(diag(Jcov(zind,zind)));

    magopt=sqrt(diag(Jocov)); %power at each sourcepoint
    [~,peakind]=max(magopt); %% index for sourcepoint with peak power


   










