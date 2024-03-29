%% SPM inversion - new neck cast - hand contraction

clear all;
close all;
clc
restoredefaultpath

mydir='D:\brainspine_data';
subjectID ='123'; %122 or %123
spmpath='D:\spm';

%whatstr='emg+abs';
whatstr='brainopt+abs';

%% paths
addpath(spmpath)
spm('defaults','EEG');
addpath(genpath('D:\brainspineconnectivity'))

%% get meta data

dataOut=getSpinalMetaData(subjectID, mydir);

filenames=dataOut.filenamesNOISE;
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
SHUFFLE=1; %this performs additional shuffled CVA at ROI
FIXORIENT=1; %empty means it calculates optimal orientation. 1 is along spinal cord
REMOVELIN=1; 
RMORTHBRAIN=contains(whatstr,'brain'); %% don't remove ortho brain if dealing with cord-muscle
FIXMODES=0;%1 means force use of all channels in source recon
CLONETRIALS=1; % use all trial data (1) rather than average (0)

freqroi=[5 35];
invtype='IID';

IMAGECROSS=0; %only image power not cross spectrum


for cnd=1%:size(filenames,1)

    % load the data
    DAll=spm_eeg_load(filenames(cnd,:));
    [a1,b1,c1]=fileparts(DAll.fullfile);
    pstr=b1(1:10); %%get left right info from filename
    cname=[a1 filesep 'clone_b1' b1 c1 ];
    D=DAll; %keep a copy

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
    
    %get index for EMG
    emgind=D.indchantype('EMG');

    %% now compute coherence and cross spectra between brain, emg, and sensors on cord
    cd(fullfile(spmpath,'external\fieldtrip\private'))

    seedperm=0;
    cohbrain=coh_meaghan(D,D.chanlabels(brainind),D.chanlabels(emgind(1)),[min(freqroi) max(freqroi)]);
    [cspect_brain,fdat_brain,trialcspect_brain]=xspectrum_meaghan(D,D.chanlabels(brainind),D.chanlabels(emgind(1)),[min(freqroi) max(freqroi)],seedperm);
    [cspect,fdat,trialcspect]=xspectrum_meaghan(D,D.chanlabels(spineind),D.chanlabels(emgind(1)),[min(freqroi) max(freqroi)],seedperm);

    fbemgind=find(contains(fdat_brain.label,D.chanlabels(emgind(1))));
    fbrainind=setdiff(1:length(fdat_brain.label),fbemgind);


%% plot power for y channels on back

yind=find(contains(fdat.grad.label,'-Y'));
chpos=fdat.grad.chanpos(yind,1:2);
chlab=fdat.grad.label(yind);
powervals=abs(fdat.fourierspctrm).^2;
%compute the power then take the mean
meanpower=squeeze(mean(powervals,1));
meanpower=mean(meanpower,2);

% in this coord frame.  , x up-down spine, y left-right, z is out of back ant-post
%spm_eeg_plotScalpData(meanpower(yind),chpos,chlab)

    %% plot some metrics of brain-emg interaction
    nchansBrain=length(cspect_brain.labelcmb); 
    
    %         figure; subplot(4,1,1);
    %         plot(cspect_brain.freq,real(cspect_brain.crsspctrm(1:nchansBrain,:))');
    %         title('real cross')
    %         subplot(4,1,2);
    %         plot(cspect_brain.freq,imag(cspect_brain.crsspctrm(1:nchansBrain,:))');
    %         title('imag cross')
    %         subplot(4,1,3)
    %         plot(cspect_brain.freq,abs(cspect_brain.crsspctrm(1:nchansBrain,:))');
    %         title('abs cross')
    %         subplot(4,1,4)
    %         plot(cohbrain.freq,cohbrain.cohspctrm(1:nchansBrain,:)');
    %         title('coh')

    %% get dominant spatial component that explains the brain-emg cross spectrum

    r1=cspect_brain.crsspctrm;
    r1=r1-mean(r1);

    [Ur,S,~]=svd(real(r1)*real(r1)'); %% use real part
    %[Ui,Si,Vi]=svd(imag(r1)*imag(r1)');
    varexp=cumsum(diag(S))./sum(diag(S));
    Nrcomp=min(find(varexp>0.99));
    Uinv=pinv(Ur); %% the spatial mixture of brain chans orthogonal to the ideal mixture for EMG-brain cross-spectrum
    replace.U=Ur(:,1); %% Ur or Uinv
    replace.ind=fbrainind;
   
    %replace one of the channels with this mixture
    if strcmp(subjectID,'136')
        replace.label='G2-35-Y';
    else
        replace.label='G2-19-Y';
    end

    %calculate coherence between ideal mixture-EMG
    [~, fdat_brainmix,~]=xspectrum_meaghan(D,D.chanlabels(brainind),D.chanlabels(emgind(1)),[min(freqroi) max(freqroi)],seedperm,[],replace);

    combs={replace.label D.chanlabels{emgind(1)}};
    cfg            = [];
    cfg.method     = 'coh';
    cfg.channelcmb = combs;

    brainmix_emg_coh      = ft_connectivityanalysis(cfg, fdat_brainmix);

    if ~isempty(layoutfile)
        load(layoutfile)

        %     cfg                  = [];
        %     cfg.parameter        = 'cohspctrm';
        %     cfg.xlim             = [15 30];
        %     cfg.refchannel       = combs{2};
        %     cfg.layout           = lay_head;
        %     cfg.interplimits='electrodes';
        %     figure; ft_topoplotER(cfg, cohbrain)
        %     colorbar
        %     title('Brain muscle coherence')


        %% plot the optimal brain mixture weights
        %to do: normalize them

        brainlabs=fdat_brain.label(~contains(fdat_brain.label,'EMG'));

        dat=spm2fieldtrip(D);
        cfg=[];
        cfg.channel=brainlabs;
        dat=ft_selectdata(cfg,dat);


        cfg=[];
        cfg.covariance='yes';
        datacov=ft_timelockanalysis(cfg,dat);

        dat=ft_timelockanalysis([],dat); %put into ft structure to plot

        dat.avg=datacov.cov*Ur(:,1);
        dat.time=dat.time(2);
        dat.var=dat.var(:,1);
        dat.dof=dat.dof(:,1);
        %% topoplot
        % cfg                  = [];
        % cfg.parameter        = 'avg';
        % cfg.layout           = lay_head;
        % cfg.interplimits='electrodes';
        % figure; ft_topoplotER(cfg, dat)
        % colorbar
        % title('opt brain mix')
    end %if layout file


    %% replace one of the chans with optimal linear mixture and compare to other channels

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

    if CLONETRIALS %put the complex data into D obj
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
    [u_src,s_src,v_src]=svd(p1'*p1); %% orientations along (longer) spine axis u_src(:,1)


    %% convert to SI units
    src_m=ft_convert_units(src,'m');
    grad_m=ft_convert_units(grad,'m');
    [subject_m] = ft_convert_units(subject, 'm');

    %% loop over forward models

    arbvals=strvcat('singlesphere'); 
    %arbvals=strvcat('singleshell','George','localspheres','singlesphere','infinite');


    allF=[]; %can loop over volume conductors and compare FE
    allR2=allF;
    allVE=allF;
    allmaxFstat=allF;

    for arbind=1:size(arbvals,1) %loop through volume conductors

        %now set up torso model

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
                %cfg.normalize           = 'yes';
                %cfg.normalizeparam     = 0.5;

                sourcemodel = ft_prepare_leadfield(cfg);
                
                src_curr=sourcemodel;
                
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

                % if we want to shift the center of the sphere----------
                 %if findstr(arbvals(arbind,:),'sphere')
                     %bmodel.o(:,3)=bmodel.o(:,3)+.1;
                     %warning('shifting sphere center')
% %                     figure;
% %                     plot3(ps(:,1),ps(:,2),ps(:,3),'c.');
% %                     hold on;
% %                     plot3(bmodel.o(:,1),bmodel.o(:,2),bmodel.o(:,3),'r*');
                 %end
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

            case 'George'
                error('This is not set up yet')
                [Lf,grad,src]=bem_leadfields_neckcast(subject,sub_fids,grad,src);
                %Lf dimensions are nchans x 3nsourcepts

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
    Jvnoise=Jvnoise./sqrt(dot(Jvnoise,Jvnoise)); %this is optimal orientation

    if ~isempty(FIXORIENT) %if user specifies orientation use this instead
        Jv=u_src(:,FIXORIENT);
    end


    % plot optimal orientation
   % plotOptOri(Jv,src,subject)
    
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
    magopt_noise=sqrt(diag(Jonoisecov));
    [~,peakind]=max(magopt); %% index for sourcepoint with peak power


    %% figure
    f1=figure;
    xvals=src.pos(:,1);
    plot(xvals,magopt,'o','MarkerFaceColor',[0 .3 .3],'MarkerEdgeColor',[0 .5 .5],'MarkerSize',5);
    ax = gca;
    ax.FontSize = 18;
    box off
    set(gcf, 'Position', [570 769 670 227]);
    xlabel('Cranial-caudal position (cm)')
    ylabel('Power (AU)')

    % Calculate the median y values for each unique x value
    unique_x = unique(xvals);
    median_y_values = splitapply(@median, magopt, findgroups(xvals));

    hold on; %
    plot(unique_x(1:end-1), median_y_values(1:end-1), 'color',[0 .5 .5], 'LineWidth', 2);

    savename=sprintf('%s %s magopt for sourceplot.pdf',subjectID,num2str(cnd));
    %exportgraphics(gcf, fullfile(figsavedir,savename), 'Resolution', 600)
   

    %plot peak power
    figure;
    plot_func_spine_dat(subject,src,magopt,grad,sens_stl_ft)
    hold on
    plot3(src.pos(peakind,1), src.pos(peakind,2),src.pos(peakind,3),'ro','MarkerSize',10,'MarkerFaceColor','r')
    savename=sprintf('%s %s %g.png',subjectID,whatstr,cnd);
    exportgraphics(gcf, fullfile(figsavedir,savename), 'Resolution', 600)

    savename=sprintf('%s %s %g_backview.png',subjectID,whatstr,cnd);

    %exportgraphics(gcf, fullfile(figsavedir,savename), 'Resolution', 600)

    %% Ybrain should be data from optimal linear mixture of channels to get emg coherence
    % Yibrain should be data from channel mixture orthogonal to this

    Ybrain=[];Yibrain=[];
    for f=1:size(fdat_brain.fourierspctrm,1)
        Ybrain(f,:)=Ur(:,1)'*squeeze(fdat_brain.fourierspctrm(f,bind,:));
        Yibrain(f,:)=Uinv(:,1)'*squeeze(fdat_brain.fourierspctrm(f,bind,:));
    end

    %% X is always spinal cord i.e. Jorient
    useJ=Jorient;
    X=Jorient;

    % Y can either be emg or brain
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

    %% linear regression: Y (EMG or brain) and ortho brain signal

    Fstatlin=zeros(length(usefreq));
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

    %% Linear regression: emg-brain

    rYbrain=zeros(size(Ybrain2));
    Fstatlin_emgbrainri=zeros(size(cind,2),1);
    pvals_emgbrainri=zeros(size(cind,2),1);


    for f=1:size(cind,2)*2 %loop through frequencies

        %Ybrain2 is copy of original brain signal, Zemg is copy of original EMG. Both are complex.

        [B,BINT,R,RINT,STATS] = regress(Ybrain2(:,f),[Zemg(:,rpcind(1,f)) Zemg(:,rpcind(2,f)) ones(size(Zemg,1),1)]);
        rYbrain(:,f)=R;
        Fstatlin_emgbrainri(f)=STATS(2);
        pvals_emgbrainri(f)=STATS(3);

    end % for f
    
   
    %save this to bootstrap variance explained later
    save(fullfile(savepath,sprintf('regdat_brainemg_%s_%s',subjectID,whatstr)), 'Ybrain2', 'Zemg','cind','rpcind')

    Fstatlin_emgbrain_real=Fstatlin_emgbrainri(cind(1,:));
    Fstatlin_emgbrain_imag=Fstatlin_emgbrainri(cind(2,:)); %% separate real and imag

    pvals_emgbrain_real=pvals_emgbrainri(cind(1,:));
    pvals_emgbrain_imag=pvals_emgbrainri(cind(2,:));

    % work out the threshold using FDR
    [~, crit_p_real, ~, adj_p_real]=fdr_bh(pvals_emgbrain_real,0.05,'dep','yes');
    [~, crit_p_imag, ~, adj_p_imag]=fdr_bh(pvals_emgbrain_imag,0.05,'dep','yes');
    realSig=find(adj_p_real<0.05);
    imagSig=find(adj_p_imag<0.05);

    fprintf('EMG brain Critical p real is %.2f and imag is %.2f\n',crit_p_real, crit_p_imag)

    %save variance removed
    linvaremgbrain_rmvd=1-var(rYbrain(:))/var(Ybrain2(:)); %residuals as prop of orig brain signal

    %plot linear results
    plotLinearMeasures(usefreq, Fstatlin_emgbrain_real, Fstatlin_emgbrain_imag, realSig, imagSig, brainmix_emg_coh,'Brain','EMG')
    savename=sprintf('lin_sub%s_emg_brain_%g.pdf',subjectID,cnd);
    exportgraphics(gcf, fullfile(figsavedir,savename), 'Resolution', 600)
   %% use residuals from linear regression to look at within frequency envelope correlations (abs)

    rYbrain_abs=abs([rYbrain(:,cind(1,:))+i*rYbrain(:,cind(2,:))]); %% put it back into complex numbers then abs
    Zemg_abs=abs([Zemg(:,cind(1,:))+i*Zemg(:,cind(2,:))]);
    
    %mean centre
    rYbrain_abs=rYbrain_abs-mean(rYbrain_abs);
    Zemg_abs=Zemg_abs-mean(Zemg_abs);

    %preallocate
    Fstat_emgbrain_env=zeros(1,numel(usefreq));
    pvals_emgbrain_env=zeros(1,numel(usefreq));
    rYbrain2=zeros([size(rYbrain_abs)]); %% save residuals

    %loop through frequencies
    for f=1:numel(usefreq)
        [B,BINT,R,RINT,STATS] = regress(rYbrain_abs(:,f),[Zemg_abs(:,f) ones(size(Zemg_abs,1),1)]);
        rYbrain2(:,f)=R;
        Fstat_emgbrain_env(f)=STATS(2);
        pvals_emgbrain_env(f)=STATS(3);
    end % for f

% work out the threshold using FDR
    [~, crit_p_env, ~, adj_p_env]=fdr_bh(pvals_emgbrain_env,0.05,'dep','yes');
    env_Sig=find(adj_p_env<0.05);

 %% plot linear results

 %put this in a function at some point!
cols=colormap(brewermap([],"Dark2"));
col3=cols(3,:); col4=cols(4,:);
col5=cols(8,:);

% Where to put significance stars
[maxValue, ~] = max(Fstat_emgbrain_env);
yValue1 = maxValue + 0.1 * maxValue; %

figure;
plot(usefreq, Fstat_emgbrain_env, 'color',col3,'LineWidth',3); hold on
if ~isempty(env_Sig)
    plot(usefreq(env_Sig), ones(1,length(env_Sig))*yValue1,'*','MarkerSize',8,'color',col3)
end
ylabel('F Statistic');
xlabel('Frequency (Hz)');
ax = gca;
ax.FontSize = 18;
ax.LineWidth=1.5; %change to the desired value     
box off
xlim([0 40])
title('EMG brain envelope correlation')


savename=sprintf('lin_ENV_sub%s_emg_brain_%g.pdf',subjectID,cnd);
%exportgraphics(gcf, fullfile(figsavedir,savename), 'Resolution', 600)

    
    %% emg-brain cva (using residuals from linear regression)
    rYbrain2=rYbrain2-mean(rYbrain2);
    CVAemgbrain_abs=spm_cva(rYbrain2,Zemg_abs);
    fprintf('EMG brain CVA chi2=%.2f p=%.3f\n',CVAemgbrain_abs.chi(1),CVAemgbrain_abs.p(1))

    if SHUFFLE
        rYbrain2=rYbrain2(randperm(size(rYbrain2,1)),:);
        rYbrain2=rYbrain2-mean(rYbrain2);
        CVA_shuf_emgbrain=spm_cva(rYbrain2,Zemg_abs);
        fprintf('CVA shuf data emg brain chi2=%.2f p=%.3f \n', CVA_shuf_emgbrain.chi(1),CVA_shuf_emgbrain.p(1))
    end


    %normalize canonical vectors
    Ncan_emgbrain=max(find(CVAemgbrain_abs.p<0.05));
    normV_emgbrain=(cov(CVAemgbrain_abs.Y))*CVAemgbrain_abs.V(:,1:Ncan_emgbrain)*inv(cov(CVAemgbrain_abs.v(:,1:Ncan_emgbrain)));
    normW_emgbrain=(cov(CVAemgbrain_abs.X))*CVAemgbrain_abs.W(:,1:Ncan_emgbrain)*inv(cov(CVAemgbrain_abs.w(:,1:Ncan_emgbrain)));

     plotNonlinearMeasures(usefreq, abs(normV_emgbrain(:,1)),abs(normW_emgbrain(:,1)),'Brain','EMG')
%     savename=sprintf('nonlin_sub%s_emg_brain_%g.pdf',subjectID,cnd);
%     exportgraphics(gcf, fullfile(figsavedir,savename), 'Resolution', 600)
%     
    %% Spinal cord and Y(EMG/brain) linear regression at source point where power is max

    Y=Y-mean(Y);
    Fstatlinri=zeros(1,numel(cind));
    pvalsri=zeros(1,numel(cind));
    rYcord=zeros([size(Y)]); %% save residuals
    X=[real(squeeze(useJ(peakind,:,:))') imag(squeeze(useJ(peakind,:,:))')]; % cord data where power is max
    X=X-mean(X);

    for f=1:size(cind,2)*2 %% want to deal with real and imag part of Y

        [B,BINT,R,RINT,STATS] = regress(Y(:,f),[X(:,rpcind(1,f)) X(:,rpcind(2,f)) ones(size(X,1),1)]);
        rYcord(:,f)=R; %% THIS WAS WRONG AS rY changed at each location on cord...

        % note also that rYcord has both real and imag components that had linear relationship removed

        Fstatlinri(f)=STATS(2);
        rsqlin(f)=STATS(1);
        pvalsri(f)=STATS(3);
    end % for f

    %save for bootstrapping variance explained in other script
    save(fullfile(savepath,sprintf('regdat_spine_%s_%s',subjectID,whatstr)), 'Y', 'X','cind','rpcind')


    Fstatlin_real=Fstatlinri(cind(1,:)); %% separate real and imag
    Fstatlin_imag=Fstatlinri(cind(2,:));

    pvals_real=pvalsri(cind(1,:)); %corresponding p values
    pvals_imag=pvalsri(cind(2,:));


    % work out the threshold using FDR
    [h, crit_p_r, ~, adj_p_real]=fdr_bh(pvals_real,0.05,'dep','yes');
    [~, crit_p_i, ~, adj_p_imag]=fdr_bh(pvals_imag,0.05,'dep','yes');
    realSig=find(adj_p_real<0.05);
    imagSig=find(adj_p_imag<0.05);

    fprintf('Critical p real is %.2f and imag is %.2f\n',crit_p_r, crit_p_i)

    %calc variance removed for point of max power
    rYMax=rYcord;
    linvar_rmvd= 1-var(rYMax(:))/var(Y(:)); %save this for point of max pwr. Y2 is a copy of Y from lin reg (is overwritten)

    %% Calculate corresponding coherence for visualization

    SC_freqdat=squeeze(useJ(peakind,:,:))'; % cord, complex
    YCoh=Y(:,cind(1,:),:) +i*Y(:,cind(2,:),:); %put y back into complex
   
    if findstr(whatstr,'brain')
        lab1='Brain';
    end

    if findstr(whatstr,'emg')
        lab1='EMG';
    end

    Coh_ft = cohxy(SC_freqdat,YCoh,fdat);


    plotLinearMeasures(usefreq, Fstatlin_real, Fstatlin_imag, realSig, imagSig, Coh_ft,lab1,'spinal cord')
    savename=sprintf('linear_sub%s_%s_cord_%g.pdf',subjectID,lab1,cnd);
    exportgraphics(gcf, fullfile(figsavedir,savename), 'Resolution', 600)
    
    %% define X Y and X0 for Y (brain or EMG)-envelope corr

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
            Ycord=abs(rYcord(:,cind(1,:),:) +i*rYcord(:,cind(2,:),:));
            X0=abs(Yibrain); % abs of orth brain : use confound here as non-linear portion may influence

        otherwise
            error('not defined')
    end % switch
%% amplitude envelope cord and emg/brain
 
%preallocate
    Fstat_env=zeros(1,numel(usefreq));
    pvals_env=zeros(1,numel(usefreq));
    rYcord2=zeros([size(Ycord)]); %% save residuals

    X=squeeze(useJ(peakind,:,:))'; % cord data

    if contains(whatstr,'abs')
        X=abs(X);
    else
        X=[real(X) imag(X)];
    end

    % here is where we need to account for/ quantify linear effect from part of cord (i.e. just before non lin stage)
    Y=squeeze(Ycord(:,:)); %% this will have had linear portion removed if required

    Y=Y-mean(Y);
    X=X-mean(X);

    %loop through frequencies
    for f=1:numel(usefreq)
        [B,BINT,R,RINT,STATS] = regress(Y(:,f),[X(:,f) ones(size(X,1),1)]);
        rYcord2(:,f)=R;
        Fstat_env(f)=STATS(2);
        pvals_env(f)=STATS(3);
    end % for f

    [~, crit_p_env, ~, adj_p_env]=fdr_bh(pvals_env,0.05,'dep','yes');
    envSig=find(adj_p_env<0.05);

% Where to put significance stars
[maxValue, ~] = max(Fstat_env);
yValue1 = maxValue + 0.1 * maxValue; %

figure;
plot(usefreq, Fstat_env, 'color',col3,'LineWidth',3); hold on
if ~isempty(envSig)
    plot(usefreq(envSig), ones(1,length(envSig))*yValue1,'*','MarkerSize',8,'color',col3)
end
ylabel('F Statistic');
xlabel('Frequency (Hz)');
ax = gca;
ax.FontSize = 18;
ax.LineWidth=1.5; %change to the desired value     
box off
xlim([0 40])
title('Cord envelope correlation')


    %% CVA Y (brain/EMG) and spinal cord, with linear part removed
    
    Y=squeeze(rYcord2(:,:)); %% use residuals from env corr
    Y=Y-mean(Y);


    X0=X0-mean(X0);
    chi2=0;
    chi2r=0;
    pval=[];
    allCVA=[];


    CVA=spm_cva(Y,X,X0);

    allCVA{1}=CVA; %for saving

    fprintf('CVA at point of peak power chi2=%.2f p=%.3f \n', CVA.chi(1),CVA.p(1))
    fprintf('Source recon VE is %.2f R2 is %.2f\n', allVE,allR2)

    Ncan=max(find(CVA.p<0.05));

    %normalize can vec
    normV=cov(CVA.Y)*CVA.V(:,1:Ncan)*inv(cov(CVA.v(:,1:Ncan)));
    normW=cov(CVA.X)*CVA.W(:,1:Ncan)*inv(cov(CVA.w(:,1:Ncan)));

    %% shuffle just for point of max pwr

    if SHUFFLE
        disp('SHUFFLING DATA !')
        Y=squeeze(rYcord2(:,:));
        rng('default') %
        Y=Y(randperm(size(Y,1)),:);
        X=squeeze(useJ(peakind,:,:))'; % cord data
        Y=Y-mean(Y);
        X=X-mean(X);
        CVA_shuf=spm_cva(Y,X,X0);
        fprintf('CVA shuf data at point of peak power chi2=%.2f p=%.3f \n', CVA_shuf.chi(1),CVA_shuf.p(1))
    end
    %% save CVA data

    cvaname=[savepath filesep sprintf('%s_%s.mat',pstr,whatstr)];
    save(cvaname,'allCVA','peakind','magopt','usefreq','CVAemgbrain_abs','Fstatlinri','Fstatlin_emgbrainri','linvar_rmvd','linvaremgbrain_rmvd');

    %% plot nonlinear features from CVA

    if contains(whatstr,'emg')
        lab1='EMG';
    else
        lab1='Brain';
    end

    plotNonlinearMeasures(usefreq, abs(normV(:,1)),abs(normW(:,1)),lab1,'Spinal Cord')
    savename=sprintf('nonlin_sub%s_%s_cord_%g.pdf',subjectID,lab1,cnd);
    exportgraphics(gcf, fullfile(figsavedir,savename), 'Resolution', 600)

    %% final CVA between precision and spinal cord
%     load(fullfile(savepath,sprintf('%s_error_all_%s',subjectID,pstr(3)))) %load trial wise errors
%     %make a matrix for CVA with ntapers repeats
%     repError = repelem(allErrors, 9);
%     X=squeeze(useJ(peakind,:,:))'; %spinal cord data
%     X=X-mean(X);
%     Y=repError';
%     Y=Y-mean(Y);
% 
% 
%     CVA_precision=spm_cva(Y,X); %tendency for 116...
% 
% 

    % to save:
    % CVA_shuf_emgbrain.p(1)
    %CVA_shuf.p(1)
    %CVA_precision.p(1)
    %Ns=size(X,1)
    %'linvar_rmvd','linvaremgbrain_rmvd'

end % for filenames



