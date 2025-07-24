%% SPM inversion - new neck cast - hand contraction

clear all;
close all;
clc
restoredefaultpath

mydir='D:\brainspine_data';
subjectID ='123'; %122 or %123
spmpath='D:\spm';
repopath='D:\brainspineconnectivity';

%whatstr='emg+abs';
whatstr='brainopt+abs';

%% paths
addpath(spmpath)
spm('defaults','EEG');
addpath(genpath(repopath))

%% get meta data

dataOut=getSpinalMetaData(subjectID, mydir);

filenames=dataOut.filenames;
posfile=dataOut.posfile;
backshape=dataOut.backshape;
cyl=dataOut.cyl;
brainchan_labels=dataOut.brainchan_labels;
badchans=dataOut.badchans;
dpath=dataOut.dpath;
savepath=dataOut.savepath;
layoutfile=dataOut.layoutfile;
sensorStl=dataOut.sensorstl;
figsavedir=dataOut.figsavedir;


if ~exist(savepath,'dir')
    mkdir(savepath)
end

if ~exist(figsavedir,'dir')
    mkdir(figsavedir)
end

cd(dpath)

%% analysis options

SHUFFLE=1; %no longer used (permutation option)
FIXORIENT=2; %which source orientation along spinal cord. Empty=>calculates optimal ori
FIXMODES=0; %1 means force use of all channels in source recon

freqroi=[1 35]; %we want to plot connectivity spectra from 1-35 Hz
freqtest=5; %but we only want to test from 5 hz
invtype='IID'; % min norm... do we want to change this?

IMAGECROSS=0; %only image power not cross spectrum


    % load the data
    DAll=spm_eeg_load(filenames(cnd,:));
    [a1,b1,c1]=fileparts(DAll.fullfile);
    %pstr=b1(1:10); %%get left right info from filename
    cname=[a1 filesep 'clone_b1' b1 c1 ];
    D=DAll; %keep a copy


    %% identify channels on spine (the ones that have positions and ori)

    grad=D.sensors('MEG');
    badgrad=badchannels(D); %bad channels
    [a1,b1,spineind]=intersect(grad.label,D.chanlabels);
    spineind=setdiff(spineind,D.badchannels); % remove bad channels
    refind=find(contains(chantype(D),'REF')); %indices for ref channesl


    %% get list of indices of channels on cortex
    brainind=[];

    for g=1:length(brainchan_labels)
        brainind=[brainind find(contains(D.chanlabels,deblank(brainchan_labels(g))))];
    end

    %% get index for EMG
    emgind=D.indchantype('EMG');

    %% compute coherence and cross spectra between brain, emg, and sensors on cord

    cd(fullfile(spmpath,'external\fieldtrip\private')) %hack because FT code is edited.
    %can probably be fixed using built in function

    seedperm=0;
    cohbrain=coh_meaghan(D,D.chanlabels(brainind),D.chanlabels(emgind(1)),[min(freqroi) max(freqroi)]);

    [cspect,fdat,trialcspect]=xspectrum_meaghan(D,D.chanlabels(spineind),D.chanlabels(emgind(1)),[min(freqroi) max(freqroi)],seedperm);
    [cspect_brain,fdat_brain,trialcspect_brain]=xspectrum_meaghan(D,D.chanlabels(brainind),D.chanlabels(emgind(1)),[min(freqroi) max(freqroi)],seedperm);

    fbemgind=find(contains(fdat_brain.label,D.chanlabels(emgind(1)))); %indices in freq dat struct for emg
    fbrainind=setdiff(1:length(fdat_brain.label),fbemgind); %indices in freq dat structure ofr brain


    %% get dominant spatial component that explains the brain-emg cross spectrum

    % we should plot this too!

    r1=cspect_brain.crsspctrm;
    r1=r1-mean(r1);

    [Ur,~,~]=svd(real(r1)*real(r1)'); %% use real part.... prob should use imag part
   
    %[Ui,Si,Vi]=svd(imag(r1)*imag(r1)');

    replace.U=Ur(:,1); %% Ur
    replace.ind=fbrainind;

    %replace one of the channels with this mixture

    if strcmp(subjectID,'136')
        replace.label='G2-35-Y';
    else
        replace.label='G2-19-Y';
    end


    %% replace one of the chans with brain component

    [cspect_brainmix,fdat_brainmix,trialcspect_brainmix]=xspectrum_meaghan(D,D.chanlabels(brainind),D.chanlabels(emgind(1)),[min(freqroi) max(freqroi)],seedperm,[],replace);
    ind=find(contains(cspect_brainmix.labelcmb,replace.label));

    %% make control signal for null space based on reference channels

    %get freq decomposition
    [~,fdat_ref,~]=xspectrum_meaghan(D,D.chanlabels(refind),D.chanlabels(emgind(1)),[min(freqroi) max(freqroi)],seedperm);

    %get rid of EMG (function above requires)
    cfg=[];
    cfg.channel=fdat_ref.label(~contains(fdat_ref.label,'EMG'));
    fdat_ref=ft_selectdata(cfg,fdat_ref);

    %mean over trials and tapers
    refdat_FC=squeeze(mean(fdat_ref.fourierspctrm,1)); %dimensions trialstapers x channels x freq

    refdat=refdat_FC-mean(refdat_FC); %mean centre

    [Uref,S,~]=svd(refdat*refdat');

    varexp=cumsum(diag(S))./sum(diag(S));

    Nrcomp=min(find(varexp>0.95));%take nr comps that explain 95% var

    refX0=[]; %multiply data by weights.

    for comps=1:Nrcomp

        for t=1:size(fdat_ref.fourierspctrm,1)

            refX0(t,comps,:)=Uref(:,comps)'*squeeze(fdat_ref.fourierspctrm(t,:,:)); %this is our reference 'control' signal

        end
    end

    %% get indices of channels within the fdat structure

    %NOTE throughout the code I think there are some redundancies with vectors
    %of indicies, could clean this up

    fspemgind=find(contains(fdat.label,D.chanlabels(emgind(1))));
    spind=setdiff(1:length(fdat.label),fspemgind);
    bind=setdiff(1:length(fdat_brain.label),fbemgind);

    femgind=find(contains(fdat.label,D.chanlabels(emgind(1))));
    fuseind=setdiff(1:length(spineind),femgind); %freqs to use

    %% average over tapers to get 1 spectrum per trial

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

%put the complex data into D obj
        
Dclone=clone(DAll,cname,[size(DAll,1) Nclonesamples size(DAll,3)]);

for f=1:Dclone.ntrials
    Dclone(spineind,cind(1,:),f)=real(squeeze(trialcspect(f,:,:)));
    Dclone(spineind,cind(2,:),f)=imag(squeeze(trialcspect(f,:,:)));
end
   
    %% load optical scan stl and sensor stls (for visualization)

    [F, V]=stlread(backshape);
    sub=[];
    sub.vertices = V;
    sub.faces = F;
    sub_red = reducepatch(sub,0.5);

    [F1, V1]=stlread(sensorStl);
    sensstl=[];
    sensstl.vertices = V1;
    sensstl.faces = F1;
    sens_stl_red = reducepatch(sensstl,0.5);

    %put this in fieldtrip structure
    subject=[];
    subject.pos=sub_red.vertices;
    subject.tri=sub_red.faces;

    sens_stl_ft=[];
    sens_stl_ft.pos=sens_stl_red.vertices;
    sens_stl_ft.tri=sens_stl_red.faces;

    %% Now set up a source space cylinder

    cyl_source=ft_read_headshape(cyl);

    res=10; %cylinder grid resolution
    src=make_spine_grid_MES_cylinder(sub_red,cyl_source,res);
    p1=src.pos;
    p1=p1-mean(p1);

    [u_src,s_src,v_src]=svd(p1'*p1); %% orientations along spine axis


    %% convert to SI units
    
    src_m=ft_convert_units(src,'m');
    grad_m=ft_convert_units(grad,'m');
    subject_m = ft_convert_units(subject, 'm');

    %% loop over forward models

    arbvals=strvcat('infinite_currentdipole');
    %arbvals=strvcat('singleshell','localspheres','singlesphere','infinite');

    allF=[]; %can loop over volume conductors and compare FE
    allR2=allF;
    allVE=allF;
    allmaxFstat=allF;

    for arbind=1:size(arbvals,1) %loop through volume conductors

        Lf=[];

        switch deblank(arbvals(arbind,:))

            case {'infinite_currentdipole'}

                cfg                     = [];
                cfg.method              = 'infinite';
                cfg.siunits=1;
                cfg.grad=grad_m;
                cfg.conductivity = 1;

                vol                     = ft_prepare_headmodel(cfg,subject_m);
                vol.type                = 'infinite_currentdipole';
                vol.unit                = 'm';

                % calculate forward
                cfg                     = [];
                cfg.sourcemodel         = src_m;
                cfg.headmodel           = vol;
                cfg.grad                = grad_m;
                cfg.reducerank          = 'no';

                sourcemodel = ft_prepare_leadfield(cfg);


                %unravel the lead fields
                startidx=1;
                for k=1:size(src.pos,1) %for each source point
                    Lf(:,startidx:startidx+2)=sourcemodel.leadfield{k};
                    startidx=startidx+3;
                end

            case {'singleshell','singlesphere','localspheres','infinite'}
                cfg=[];
                cfg.grad=grad_m;
                cfg.conductivity = 1;
                cfg.units='m';
                cfg.siunits=1;
                cfg.method=deblank(arbvals(arbind,:));
                [bmodel] = ft_prepare_headmodel(cfg,subject_m);

                bmodel.unit='m';

                cfg = [];
                cfg.headmodel = bmodel;
                cfg.sourcemodel = src_m;
                cfg.reducerank= 'no';
                cfg.grad=grad_m;
                cfg.siunits=1;

                [sourcemodel] = ft_prepare_leadfield(cfg);


                %unravel leadfields
                startidx=1;
                for k=1:size(src.pos,1) %for each source point
                    Lf(:,startidx:startidx+2)=sourcemodel.leadfield{k};
                    startidx=startidx+3;
                end


            otherwise
                error('No forward model defined')

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

        Dout=spm_eeg_invert_classic_volumetric(D,1);

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
    end % for arbind (loop through volume conductors)

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
    Jnoise=zeros(Nx*3,length(fcind),ntapers);

    Jcov=zeros(Nx*3,Nx*3);
    Jnoisecov=zeros(Nx*3,Nx*3);

    for f=1:ntapers
        Jtrial=M*squeeze(spdata(f,:,fcind));
        Jnoisetr=M*eye(size(squeeze(spdata(f,:,fcind))));
        Jcov=Jcov+cov(Jtrial');
        J(:,:,f)=Jtrial;
        Jnoise(:,:,f)=Jnoisetr;
        Jnoisecov=Jnoisecov+cov(Jnoisetr');
    end
    Jcov=Jcov./ntapers;
    Jnoisecov=Jnoisecov./ntapers;

    xind=1:3:Nx*3; %indices for x y and z oriented sources
    yind=2:3:Nx*3;
    zind=3:3:Nx*3;

    %% get optimal orientation based on whole cord

    Jv=sqrt([sum(diag(Jcov(xind,xind))) sum(diag(Jcov(yind,yind))) sum(diag(Jcov(zind,zind)))]);
    Jv=Jv./sqrt(dot(Jv,Jv)); %this is optimal orientation

    %get orientation of noise
    Jvnoise=sqrt([sum(diag(Jnoisecov(xind,xind))) sum(diag(Jnoisecov(yind,yind))) sum(diag(Jnoisecov(zind,zind)))]);
    Jvnoise=Jvnoise./sqrt(dot(Jvnoise,Jvnoise));

    if ~isempty(FIXORIENT) %if user specifies orientation use this instead
        Jv=u_src(:,FIXORIENT);
    end


    % plot orientation of sources
    %     figure
    %     plotOptOri(Jv,src,subject)
    %     title('Source orientation')

    Jorient=J(xind,:,:).*Jv(1)+J(yind,:,:).*Jv(2)+J(zind,:,:).*Jv(3); %source data for opt orientation or fixorient
    Jnoiseorient=Jnoise(xind,:,:).*Jv(1)+Jnoise(yind,:,:).*Jv(2)+Jnoise(zind,:,:).*Jv(3);
    Jocov=zeros(length(xind),length(xind));
    Jonoisecov=zeros(length(xind),length(xind));

    for f=1:ntapers
        Jtrial=Jorient(:,:,f);
        Jnoisetrial=Jnoiseorient(:,:,f);
        Jocov=Jocov+cov(Jtrial');
        Jonoisecov=Jonoisecov+cov(Jnoisetrial');
    end

    Jocov=Jocov./ntapers;
    Jonoisecov=Jonoisecov./ntapers;

    magopt=sqrt(diag(Jocov)); %power at each sourcepoint
    magopt_noise=sqrt(diag(Jonoisecov)); %same for noise
    [~,peakind]=max(magopt); %% index for sourcepoint with peak power

    %% Make source figure
    f1=figure;
    xvals=src.pos(:,1);
    plot(xvals,magopt,'o','MarkerFaceColor',[0 .3 .3],'MarkerEdgeColor',[0 .5 .5],'MarkerSize',5);
    ax = gca;
    ax.FontSize = 18;
    ax.LineWidth=1.25;
    box off
    set(gcf, 'Position', [570 769 670 227]);

    cmin = min(magopt(:));
    cmax = max(magopt(:));

    % Calculate the median y values for each unique x value
    unique_x=unique(xvals);
    median_y_values=splitapply(@median, magopt, findgroups(xvals));

    hold on;
    plot(unique_x(1:end-1), median_y_values(1:end-1), 'color',[0 .5 .5], 'LineWidth', 2);

    savename=sprintf('%s %s magopt for sourceplot.pdf',subjectID,num2str(cnd));
    exportgraphics(gcf, fullfile(figsavedir,savename), 'Resolution', 600)

    %% plot peak power

    figure;
    plot_func_spine_dat(subject,src,magopt,grad,sens_stl_ft)
    hold on

    if strcmp(subjectID,'116') && cnd==1
        warning('Hard coded to find local max for 116')
        plot3(src.pos(peakind,1), src.pos(peakind,2),src.pos(peakind,3),'ro','MarkerSize',10)

        [~,I]=maxk(magopt,6);

        peakind=I(6);
        plot3(src.pos(peakind,1), src.pos(peakind,2),src.pos(peakind,3),'ro','MarkerSize',10,'MarkerFaceColor','r')


    elseif strcmp(subjectID,'116') && cnd==2
        plot3(src.pos(peakind,1), src.pos(peakind,2),src.pos(peakind,3),'ro','MarkerSize',10)
        warning('Hard coded to find local max for 116')

        [~,I]=maxk(magopt,7);
        peakind=I(7);
        plot3(src.pos(peakind,1), src.pos(peakind,2),src.pos(peakind,3),'ro','MarkerSize',10,'MarkerFaceColor','r')

    end


    plot3(src.pos(peakind,1), src.pos(peakind,2),src.pos(peakind,3),'ro','MarkerSize',10,'MarkerFaceColor','r')
    savename=sprintf('%s %s %g.png',subjectID,whatstr,cnd);
    exportgraphics(gcf, fullfile(figsavedir,savename), 'Resolution', 600)

    %    savename=sprintf('%s %s %g_backview.png',subjectID,whatstr,cnd);
    %    exportgraphics(gcf, fullfile(figsavedir,savename), 'Resolution', 600)

    %% Ybrain should be data from optimal linear mixture of channels to get cross spectrum with emg

    Ybrain=[];
    for f=1:size(fdat_brain.fourierspctrm,1)
        Ybrain(f,:)=Ur(:,1)'*squeeze(fdat_brain.fourierspctrm(f,bind,:));
    end

    %% EMG-Brain Functional connectivity

    %mean centre
    Ybrain_complex=Ybrain-mean(Ybrain);
    emgdat_complex=emdata-mean(emdata);

    freqtestidx=find(usefreq==freqtest); %we only want to test functional connectivity from % hz
    Prec=[]; %do not include precision here

    whichanalysis='emgbrain';
    [chi_lin_emgbr, chi_env_emgbr, chi_cross_emgbr,df_lin_emgbr, df_env_emgbr, df_cross_emgbr] ...
        = new_CVA_FC2(Ybrain_complex,emgdat_complex,refX0,Prec,usefreq,freqtestidx,whichanalysis,figsavedir,subjectID,cnd);

    %compute coherence for comparison
    %XY_coh= cohxy(emgdat_complex,Ybrain_complex,fdat);

    %     figure
    %     plot(XY_coh.freq, squeeze(XY_coh.cohspctrm(1,2,:)))
    %% Spinal cord and Y(EMG/brain) functional connectivity

    % Y can either be emg or brain
    Y=[];

    if findstr(whatstr,'brain')
        Y=Ybrain;
        whichanalysis='braincord';
    end

    if findstr(whatstr,'emg')
        Y=emdata;
        whichanalysis='cordemg';
    end

    if isempty(Y)
        error('Y not defined');
    end


    % X is always spinal cord i.e. Jorient
    useJ = Jorient; %spinal cord data at all source points
    X=squeeze(useJ(peakind,:,:))'; % cord data where power is max
    X=X-mean(X);
    Prec=[];
    [chi_lin_spine, chi_env_spine, chi_cross_spine, df_lin_spine, df_env_spine, df_cross_spine] ...
        = new_CVA_FC2(Y,X,refX0,Prec,usefreq,freqtestidx,whichanalysis,figsavedir,subjectID,cnd);

    % XY_coh = cohxy(X,Y,fdat);

    %figure
    %plot(XY_coh.freq, squeeze(XY_coh.cohspctrm(1,2,:)))

    %% now add precision as a column to design matrix and run cvas again

    %load Error for all one second periods
    load(fullfile(savepath,sprintf('%s_error_all_%s',subjectID,pstr(3)))) %load trial wise errors

    tapspertrial=ntapers/size(D,3); %repeat for number of tapers
    errorsrep = repelem(allErrors, tapspertrial);

    %mean centre
    errorsrep=errorsrep-mean(errorsrep);
    Prec=errorsrep';


    %spinal cord and brain or EMG WITH PRECISION
    [chi_lin_P_spine, chi_env_P_spine, chi_cross_P_spine,df_lin_P_spine, df_env_P_spine, df_cross_P_spine] ...
        =new_CVA_FC2(Y,X,refX0,Prec,usefreq,freqtestidx,whichanalysis,figsavedir,subjectID,cnd);


    %EMG brain WITH PRECISION
    whichanalysis='emgbrain';
    [chi_lin_P_emgbr, chi_env_P_emgbr, chi_cross_P_emgbr, df_lin_P_emgbr, df_env_P_emgbr, df_cross_P_emgbr] ...
        = new_CVA_FC2(Ybrain_complex,emgdat_complex,refX0,Prec,usefreq,freqtestidx,whichanalysis,figsavedir,subjectID,cnd);

    % ratio of chi square values norm by df is F stat
    Fvals_lin_emgbr=(chi_lin_P_emgbr./df_lin_P_emgbr)./(chi_lin_emgbr./df_lin_emgbr);
    Fvals_env_emgbr=(chi_env_P_emgbr./df_env_P_emgbr)./(chi_env_emgbr./df_env_emgbr);
    Fvals_cross_emgbr=(chi_cross_P_emgbr/df_cross_P_emgbr)/(chi_cross_emgbr/df_cross_emgbr);

    Fvals_lin_spine=(chi_lin_P_spine./df_lin_P_spine)./(chi_lin_spine./df_lin_spine);
    Fvals_env_spine=(chi_env_P_spine./df_env_P_spine)./(chi_env_spine./df_env_spine);
    Fvals_cross_spine=(chi_cross_P_spine/df_cross_P_spine)/(chi_cross_spine/df_cross_spine);

    
    %% calculate p values from F stats and df and correct for mult comp
    %could put this in a function to neaten up
    
    p_values_lin_emgbr = 1 - fcdf(Fvals_lin_emgbr, df_lin_P_emgbr, df_lin_emgbr);
    
    [~, ~, ~, adj_p_lin_emgbr]=fdr_bh(p_values_lin_emgbr(freqtestidx:end),0.05,'dep','no');
    
    lin_emgbrain_Sig=find(adj_p_lin_emgbr<0.05);


    
    if ~isempty(lin_emgbrain_Sig)
        fprintf('sig change, check out which freqs\n')
    else
        fprintf('no sig change for lin emg brain\n')
    end


    if sum(p_values_lin_emgbr<0.05)> 0
        idx=find(p_values_lin_emgbr<0.05);
        fprintf('P value %g for lin emg brain sig without MC correction at:\n',p_values_lin_emgbr(idx))
        disp(usefreq(idx))
    end

    %% same for env

    p_values_env_emgbr = 1 - fcdf(Fvals_env_emgbr, df_env_P_emgbr, df_env_emgbr);
   
    [~, ~, ~, adj_p_env_emgbr]=fdr_bh(p_values_env_emgbr(freqtestidx:end),0.05,'dep','no');
    
    env_emgbrain_Sig=find(adj_p_env_emgbr<0.05);

    if ~isempty(env_emgbrain_Sig)
        fprintf('sig change, check out which freqs\n')
    else
        fprintf('no sig change for env emg brain\n')
    end

    if sum(p_values_env_emgbr<0.05)> 0
        idx=find(p_values_env_emgbr<0.05); %indexed from test freq
        fprintf('P value %g for env emg brain sig without MC correction at:\n',p_values_env_emgbr(idx))
        disp(usefreq(idx))
    end

   
    %% one p value for cross freq cva comparison

    p_value_cross_emgbr = 1 - fcdf(Fvals_cross_emgbr, df_cross_P_emgbr, df_cross_emgbr);
    
    fprintf('p value for emgbr cross freq is %g\n',p_value_cross_emgbr)
   
    %% now do the same for spinal cord CVA comparisons
    %phase
    p_values_lin_spine = 1 - fcdf(Fvals_lin_spine, df_lin_P_spine, df_lin_spine);
    
    [~, crit_p_lin_spine, ~, adj_p_lin_spine]=fdr_bh(p_values_lin_spine(freqtestidx:end),0.05,'dep','no');
    
    lin_spine_Sig=find(adj_p_lin_spine<0.05);

     if ~isempty(lin_spine_Sig)
        fprintf('sig change, check out which freqs\n')
    else
        fprintf('no sig change for lin spine\n')
     end

    if sum(p_values_lin_spine<0.05)> 0
        idx=find(p_values_lin_spine<0.05);
        fprintf('P value %g for lin spine sig without MC correction at:\n',p_values_lin_spine(idx))
        disp(usefreq(idx))
    end

    %envelope
    p_values_env_spine = 1 - fcdf(Fvals_env_spine, df_env_P_spine, df_env_spine);
    
    [~, ~, ~, adj_p_env_spine]=fdr_bh(p_values_env_spine(freqtestidx:end),0.05,'dep','no');
    
    env_spine_Sig=find(adj_p_env_spine<0.05);


   % pvalues_fdr = mafdr(p_values_env_spine, 'BHFDR', true);


     if ~isempty(env_spine_Sig)
        fprintf('sig change, check out which freqs\n')
    else
        fprintf('no sig change for env spine\n')
     end

    if sum(p_values_env_spine<0.05)> 0
        idx=find(p_values_env_spine<0.05);
        fprintf('P value %g for env spine sig without MC correction at:\n',p_values_env_spine(idx))
        disp(usefreq(idx))
    end


    p_value_cross_spine = 1 - fcdf(Fvals_cross_spine, df_cross_P_spine, df_cross_spine);
   
    fprintf('p value for spine cross freq is %g\n',p_value_cross_emgbr)





