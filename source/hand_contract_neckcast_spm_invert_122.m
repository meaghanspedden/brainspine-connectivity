%% SPM inversion - new neck cast - hand contraction July 2023

clear all;
close all;

%filenames for concatenated runs for right hand and left hand
filenames=strvcat('D:\MSST001\sub-OP00122\ses-001\meg\pprhandoe1000msfdfflo45hi5ds_sub-OP00122_ses-001_task-static_right_run-002_meg.mat',...
        'D:\MSST001\sub-OP00122\ses-001\meg\pplhandoe1000msfdfflo45hi5ds_sub-OP00122_ses-001_task-static_left_run-001_meg.mat');

posfile='D:\MSST001\sub-OP00122\ses-001\meg\ds_sub-OP00122_ses-001_positions.tsv';
backshape='D:\OP00122_experiment\cast space\OP0015_seated.stl'; %stl in same space as sensors
cyl='D:\OP00122_experiment\cast space\cylinder_good_space.stl'; %cylinder for source space


dpath='D:\MSST001\sub-OP00122\ses-001\meg'; %data path
savepath='D:\MSST001\Coh_results00122'; %save path

if ~exist(savepath,'dir')
    mkdir(savepath)
end


brainchan_labels={'19','DG', 'OH', 'A1','1B', 'A9','JS','A6','DJ','MY','DS','OK','MI','17'};
badchans={'ML-X','ML-Y','ML-Z','K4-Z','K4-X'}; %% 

%%
addpath D:\torso_tools
addpath D:\analyse_OPMEG
addpath D:\spm12

spm('defaults','EEG');
addpath(genpath('D:\brainspineconnectivity'))


%% analysis options

SHUFFLE=0; %for permutation test
allcanfilenames=[];


whatstr='emg+abs';
%whatstr='brainopt+abs';
%whatstr='orthbrain+brainopt';
%whatstr='emg';

FIXORIENT=[];%[1 0 0]; % fix orientation pointing up
REMOVELIN=1; %remove linear part from CVA
RMORTHBRAIN=contains(whatstr,'brain'); %% remove orthogonal brain mixture
FIXMODES=0; %force use of all channels in source recon
CLONETRIALS=1; % use all trial data (1) rather than average (0)


freqroi=[5 30];

%invtype='EBBr95';
invtype='IID';
%invtype='EBB';

IMAGECROSS=0;
COMB3AXES=1;

sub_fids=[];


%For back scan 116
% sub_fids= [-201 238 123; %r shoulder
%     -196 -207 157; %l shoulder
%     -221 8.92 219.038  ]; %sternum


%% for back scan 115, corresponds to experiment 122

sub_fids= [19.131  -1286.201  241.474; %r shoulder
          24.208  -1750.632 262.910; %l shoulder
         -63.22  -1535.591 371.745]; %sternum




for cnd=2 %:size(filenames,1),

    DAll=spm_eeg_load(filenames(cnd,:));
    [a1,b1,c1]=fileparts(DAll.fullfile);
    pstr=b1(1:10); %%get l/r info from filename

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
    cd  D:\fieldtrip-master\fieldtrip-master\private

    seedperm=0;
    cohbrain=coh_meaghan(D,D.chanlabels(brainind),D.chanlabels(emgind(1)),[min(freqroi) max(freqroi)]);
    [cspect_brain,fdat_brain,trialcspect_brain]=xspectrum_meaghan(D,D.chanlabels(brainind),D.chanlabels(emgind(1)),[min(freqroi) max(freqroi)],seedperm);
    [cspect,fdat,trialcspect]=xspectrum_meaghan(D,D.chanlabels(spineind),D.chanlabels(emgind(1)),[min(freqroi) max(freqroi)],seedperm);

    fbemgind=find(contains(fdat_brain.label,D.chanlabels(emgind(1))));
    fbrainind=setdiff(1:length(fdat_brain.label),fbemgind);

    %% plot some metrics of brain-emg interaction
    nchansBrain=length(cspect_brain.labelcmb); %

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
    [Ur,S,V]=svd(real(r1)*real(r1)');
    [Ui,Si,Vi]=svd(imag(r1)*imag(r1)');

    Uinv=pinv(Ur); %% the spatial mixture of brain chans orthogonal to the ideal mixture for EMG-brain cross-spectrum
   
    %CHECKMIX=1;
    
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
    spind=setdiff(1:length(fdat.label),fspemgind); %spinal cord
    bind=setdiff(1:length(fdat_brain.label),fbemgind); %brain

    femgind=find(contains(fdat.label,D.chanlabels(emgind(1)))); %emg
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
        trialcspect=trialfspect; %fourier spectrum
    end

    Nclonesamples=length(cspect.freq)*2; % x2 for real and imag part
    cind=[];
    cind(1,:)=1:(Nclonesamples/2);
    cind(2,:)=((Nclonesamples/2)+1):Nclonesamples;

    imagefreq=cspect.freq;

    if CLONETRIALS %use all trials, real and imag
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
    %% read and reduce torso stl

    [F, V]=stlread(backshape);
    sub=[];
    sub.vertices = V;
    sub.faces = F;
    sub_red = reducepatch(sub,0.5);

    %put this in fieldtrip structure
    subject=[];
    subject.pos=sub_red.vertices;
    subject.tri=sub_red.faces;


    %% source reconstruction

    %save metrics for different forward models
    allF=[];
    allR2=allF;
    allVE=allF;
    allmaxFstat=allF;

    %option to loop through different forward models
    %arbvals=strvcat('singleshell','George','localspheres','singlesphere','infinite');
    
    arbvals=strvcat('infinite');

    cyl_source=ft_read_headshape(cyl);


        for arbind=1:size(arbvals,1)

            count=0;

            res=10; %resolutin for source grid in cylinder

            [src]=make_spine_grid_MES_cylinder(sub_red,cyl_source,res);
            

            %% now set up torso model
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
                    [sourcemodel] = ft_prepare_leadfield(cfg, grad); 
                    
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
                D.inv{1}.inverse.A=[]; 
            end
            
            D.inv{1}.inverse.no_temporal_filter=1;
            D.inv{1}.inverse.complexind=cind;
            D.inv{1}.forward.modality='MEG';
            D.inv{1}.forward.sensors=grad;
            D.inv{1}.forward.siunits=1; %

            D.val=1;
            D.inv{1}.inverse.Nt=[] ;


            Dout=spm_eeg_invert_classic_volumetric(D,1,Lf); %volumetric MSP

            %% mapping between source and sensor space M

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

%plotOptOri([1 0 0],src,subject)


if ~isempty(FIXORIENT) %overwrite optimal orientation with fixorient if not empty
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

% [maxval maxidx]=max(magopt);
% srctest=src;
% magopt(maxidx)=[];
% srctest.pos(maxidx,:)=[]; srctest.inside(maxidx)=[];
% 
figure
plot_func_spine_dat(subject,src,magx,grad)
title('X oriented sources')
% % 
figure
plot_func_spine_dat(subject,src,magz,grad)
title('Z oriented sources')
% 
figure
plot_func_spine_dat(subject,src,magy,grad)
title('Y oriented sources')




J=J-mean(J,3);
%Jorient=Jorient-mean(Jorient,3); %% along one axis

Ybrain=[];Yibrain=[];
for f=1:size(fdat_brain.fourierspctrm,1)
    Ybrain(f,:)=Ur(:,1)'*squeeze(fdat_brain.fourierspctrm(f,bind,:)); %spatail mixture explaining brain-muscle cross spectrum
    Yibrain(f,:)=Uinv(:,1)'*squeeze(fdat_brain.fourierspctrm(f,bind,:)); %orthogonal mixture
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

Ybrain2=[real(Ybrain) imag(Ybrain)];

Zemg=[real(emdata) imag(emdata)]; %% may need this later

if findstr(whatstr,'emg')
    Y=[real(emdata) imag(emdata)];
end

if isempty(Y)
    error('Y not defined');
end

%% now optionally remove linear parts
Fstatlin=zeros(length(usefreq),maxk);
Yconfound=[real(Yibrain) imag(Yibrain)];

%% removing orthogonal brain signal
Yr1=zeros(size(Y));
rpcind=[cind cind];
for f=1:size(cind,2)*2,
    [B,BINT,Yr1(:,f),RINT,STATS] = regress(Y(:,f),[Yconfound(:,rpcind(1,f)) Yconfound(:,rpcind(2,f)) ones(size(Yconfound,1),1)]);
    Forthri(f)=STATS(2);
end;
Forth=(Forthri(cind(1,:))+Forthri(cind(2,:)))/2; %% combine real and imag with average

if RMORTHBRAIN,
    Y=Yr1;
end;
%% get linear emg-brain
rYbrain=zeros(size(Ybrain2));

for f=1:size(cind,2)*2,
    [B,BINT,R,RINT,STATS] = regress(Ybrain2(:,f),[Zemg(:,rpcind(1,f)) Zemg(:,rpcind(2,f)) ones(size(Zemg,1),1)]);
    rYbrain(:,f)=R;

    Fstatlin_emgbrainri(f)=STATS(2);

end; % for f
% now non-linear brain emg
Fstatlin_emgbrain=(Fstatlin_emgbrainri(cind(1,:))+Fstatlin_emgbrainri(cind(2,:)))/2; %% combine real and imag with average
linvaremgbrain_rmvd=1-var(rYbrain(:))/var(Ybrain2(:));
rYbrain=abs([rYbrain(:,cind(1,:))+i*rYbrain(:,cind(2,:))]); %% put it back into complex numbers then abs
Zemg=abs([Zemg(:,cind(1,:))+i*Zemg(:,cind(2,:))]);
Ybrain=Ybrain-mean(Ybrain);
Zemg=Zemg-mean(Zemg);
CVAemgbrain=spm_cva(rYbrain,Zemg);


rY=zeros(size(Y)); % brain or emg
rX=rY;
Y=Y-mean(Y);
for k=1:maxk,
    k
    X=[real(squeeze(useJ(k,:,:))') imag(squeeze(useJ(k,:,:))')]; % cord


    X=X-mean(X);

    for f=1:size(cind,2)*2, %% want to deal with real and imag part of Y

        [B,BINT,R,RINT,STATS] = regress(Y(:,f),[X(:,rpcind(1,f)) X(:,rpcind(2,f)) ones(size(X,1),1)]);
        rY(:,f)=R;
        Fstatlinri(f,k)=STATS(2);
        rsqlin(f,k)=STATS(1);
    end; % for f

end; % for k

Fstatlin=(Fstatlinri(cind(1,:),:)+Fstatlinri(cind(2,:),:))/2; %% combine real and imag with average
if REMOVELIN,
    warning('removing linear predictions over freq')
    linvar_rmvd=1-var(rY(:))/var(Y(:));
    fprintf('\n removing %3.2f percent variance',linvar_rmvd*100);
    Y=rY;

end;



switch whatstr %% flags  RMORTHBRAIN and REMOVELIN now take care of options

    case 'brainopt',

        X0=[]; %% already removed by RMORTH flag
    case 'brainopt+abs',

        Y=abs([Y(:,cind(1,:))+i*Y(:,cind(2,:))]); %% put it back into complex numbers then abs
        X0=abs(Yibrain); % % abs of orth brain : use confound here as non-linear portion may influence

    case 'orthbrain+brainopt',
        if RMORTHBRAIN,
            error('orthbrain already removed as confound');
        end;
        X0=Y; %% null space is optimal brain signal
        Y=[real(Yibrain) imag(Yibrain)]; %% signal space is orthogonal (to optimal) brain signal

    case 'emg',

        X0=[];
    case 'emg+abs'
        Yc=abs([Y(:,cind(1,:))+i*Y(:,cind(2,:))]); %% put it back into complex numbers then abs
        Y=abs(Yc);

        X0=abs(Yibrain); % abs of orth brain : use confound here as non-linear portion may influence

    otherwise
        error('not defined')
end; % switch

if size(X)~=size(Y),
    error('X and Y don not match')
end;

if ~isempty(X0)
    if size(X0)~=size(X),
        error('X and X0 do not match')
    end;
end;

Y=Y-mean(Y);
if SHUFFLE,
    warning('SHUFFLING DATA !')
    Y=Y(randperm(size(Y,1)),:);
end;

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




    X=X-mean(X);

    CVA=spm_cva(X,Y,X0);
    chi2(k)=CVA.chi(1);
    pval(k)=CVA.p(1);
    cvar2(k)=CVA.r(1).^2;
    allCVA{k}=CVA;
end; % for k


if COMB3AXES


   [~,peakind]=max(magopt); %% peak power


    figure;
    plot_func_spine_dat(subject,src,magopt,grad)
    title('power')
    hold on
    plot3(src.pos(peakind,1), src.pos(peakind,2),src.pos(peakind,3),'ro','MarkerSize',12)


    [dum,peakindchi]=max(chi2); %% peak stat
    figure;
    plot_func_spine_dat(subject,src,chi2,grad)
    title('chi2')
    hold on
    plot3(src.pos(peakindchi,1), src.pos(peakindchi,2),src.pos(peakindchi,3),'ro','MarkerSize',12)
    title('CVA chi2')

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
    title('Brain','FontSize',20)
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
            h=plot(usefreq,Fstatlin(1:length(usefreq),peakind),usefreq,Forth);
            set(h(1),'Linewidth',4);
            if RMORTHBRAIN
                legend('Linear post removal','Brain orth removed');
            else
                legend('Linear','Brain orth not removed');
            end;
    
            ylabel('Fstat')
            xlabel('freq')
    %
    %         title(sprintf('Linear interactions %s ',pstr))




    figure;
    plot(usefreq,Fstatlin(1:length(usefreq),peakind),'LineWidth',3,'color','k');

    ylabel('Fstat','Fontsize',20)
    xlabel('Frequency (Hz)','FontSize',20)
    ax = gca;
    ax.FontSize = 20;
    %legend('Left contraction','Right contraction','FontSize',20)
    box off
    %ylim([0 2.5])

    savename=['D:\figsfortalk\linear_right_brainspine.png'];

    %exportgraphics(gcf, savename, 'Resolution', 600)

    %%
    if Ncan>0,
        plotVw=zeros(maxk,length(normW));
        plotWw=plotVw;
        for k=1:maxk,
            CVA=allCVA{k};
            rsquare=CVA.cva(1:Ncan);
            plotV=abs((cov(CVA.Y))*CVA.V(:,1:Ncan));
            plotW=abs((cov(CVA.X))*CVA.W(:,1:Ncan));
            plotVw(k,:)=plotV*rsquare; %% weighted by variance explained
            plotWw(k,:)=plotW*rsquare;
        end;
        %             figure;
        %             subplot(1,2,1);
        %             h=imagesc(usefreq,lzrange,plotVw)
        %             set(gca,'Ydir','normal')
        %             title(sprintf('Non-Linear interactions: '))
        %             ylabel('height');
        %             xlabel('freq')
        %
        %             subplot(1,2,2);
        %             h=imagesc(usefreq,lzrange,plotWw)
        %             set(gca,'Ydir','normal')
        %             title(sprintf('%s to cord %s',whatstr,pstr))
        %             ylabel('height');
        %             xlabel('freq')


    end % if Ncan>0
    
    subplot(1,2,1);
    imagesc(usefreq,lxrange,plotVw)
    set(gca,'Ydir','normal')
    title('Component 1','FontSize',12)
    %title(sprintf('Non-Linear interactions: '))
    ylabel('Height up spinal cord (mm)','FontSize',20);
    xlabel('Frequency (Hz)','FontSize',20)
    ax = gca;
    ax.FontSize = 20;
    %colormap(brewermap([],"OrRd"))


    subplot(1,2,2);
    h=imagesc(usefreq,lxrange,plotWw)
    set(gca,'Ydir','normal')
    title('Component 2','FontSize',12)

    %title('W')
    %title(sprintf('%s to cord %s',whatstr,pstr))
    ylabel('Height up spinal cord (mm)','FontSize',20);
    xlabel('Frequency (Hz)','FontSize',20)
    ax = gca;
    ax.FontSize = 20;

    %savename=['D:\figsfortalk\nonlinrightemg.png'];

    %exportgraphics(gcf, savename, 'Resolution', 600)

else % ~COMB3AXES
    %
    %         subplot(3,1,1)
    %         plot(lzrange,chi2(xind));
    %         ylabel('left-right')
    %         legend(num2str([min(pval(xind))]));
    %         title(sprintf('%s,%s',pstr,whatstr))
    %
    %         subplot(3,1,2)
    %         plot(lzrange,chi2(zind));
    %         ylabel('up-down')
    %         legend(num2str([min(pval(zind))]));
    %         subplot(3,1,3)
    %         plot(lzrange,chi2(yind));
    %         ylabel('in-out')
    %         xlabel('height up cord')
    %         legend(num2str([min(pval(yind))]));

    %% now look at which frequency components at that height
    figure;
    for orientmode=1:3,
        switch orientmode,
            case 1,
                lookind=xind;
                ostr='l-r';
            case 2,
                lookind=zind;
                ostr='u-d'
            case 3
                lookind=yind;
                ostr='in-out';
        end;
        chiplotW=[];
        chiplotV=[];
        for g=1:length(lookind),
            g
            CVA=allCVA{lookind(g)};
            normV=cov(CVA.Y)*CVA.V;
            normW=cov(CVA.X)*CVA.W;
            chival=CVA.chi(1);
            Ncan=max(find(CVA.p<0.05));

            %chiplotW(g,:)=abs(sqrt(chival)*normW(:,1)./max(normW(:,1)));
            %chiplotW(g,:)=sqrt(chival)*sum(abs(normW(:,1:Ncan)),2);
            chiplotW(g,:)=sum(abs(normW(:,1:Ncan)),2);
            %chiplotV(g,:)=abs(sqrt(chival)*normV(:,1)./max(normV(:,1)));
            %chiplotV(g,:)=sqrt(chival)*sum(abs(normV(:,1:Ncan)),2);
            chiplotV(g,:)=sum(abs(normV(:,1:Ncan)),2);

        end;
    end

    %             subplot(3,3,orientmode*3-2);
    %
    %             imagesc(freqroi,lzrange,chiplotW);
    %             set(gca,'Ydir','normal')
    %             title([ostr ' W'])
    %             colorbar;
    %             subplot(3,3,orientmode*3-1);
    %
    %             imagesc(freqroi,lzrange,chiplotV);
    %             set(gca,'Ydir','normal')
    %             title([ostr ' V'])
    %             colorbar;
    %             subplot(3,3,orientmode*3);
    %
    %             imagesc(freqroi,lzrange,chiplotV.*chiplotW);
    %             set(gca,'Ydir','normal')
    %             title([ostr ' W*V'])
    %             colorbar;
    %
    %         end; % for orient mode

    figure;
    %%



    %% now look at which frequency components at that height
    for orientmode=1:3,
        switch orientmode,
            case 1,
                lookind=xind;
                ostr='l-r';
            case 2,
                lookind=zind;
                ostr='u-d'
            case 3
                lookind=yind;
                ostr='in-out';
        end;

        [val,maxind]=max(chi2(lookind))
        chi2(yind(maxind))
        CVA=allCVA{lookind(maxind)};
        Ncan=max(find(CVA.p<0.05)); %% take significant canonical modes
        canvarX=CVA.w(:,1:Ncan);
        canvarY=CVA.v(:,1:Ncan);
        canfilename=[D.path filesep 'canvar_' ostr '_' whatstr '_' D.fname];
        save(canfilename,'CVA','canvarX','canvarY');
        allcanfilenames=strvcat(allcanfilenames,canfilename);
        chi2roi=CVA.chi(1);
        pvalroi=CVA.p(1);

        %% SO CVA.W(:,1) etc are the canonical vectors for the X
        % CVA.w is the canonical variate (=X*CVA.W)

        %% SO CVA.V(:,1) etc are the canonical vectors for the Y
        % CVA.v is the canonical variate (=Y*CVA.V)

        normV=cov(CVA.Y)*CVA.V;
        normW=cov(CVA.X)*CVA.W;


        if Ncan>0,
            realind=1:length(usefreq);
            imind=length(usefreq)+1:2*length(usefreq);

            figure;
            if ~contains(whatstr,'abs') && ~contains(whatstr,'angle'),
                subplot(2,1,1)
                h=plot(usefreq,sum(abs(normV(realind,1:Ncan)),2),'g',usefreq,sum(abs(normV(imind,1:Ncan)),2),'r');
                title(sprintf('Height %d, Sum %d modes, %s, %s',round(lzrange(maxind)),Ncan,whatstr,ostr))
                subplot(2,1,2)
                h=plot(usefreq,sum(abs(normW(realind,1:Ncan)),2),'g',usefreq,sum(abs(normW(imind,1:Ncan)),2),'r');
                title(sprintf('Sum %d modes, cord',Ncan))
            else % abs
                subplot(2,1,1)
                h=plot(usefreq,sum(abs(normV(realind,1:Ncan)),2),'m');
                title(sprintf('Height %d, Sum %d modes, %s, %s',round(lzrange(maxind)),Ncan,whatstr,ostr))
                hold on;

                h2=plot(usefreq,abs(normV(realind,1:Ncan)));


                subplot(2,1,2)
                h=plot(usefreq,sum(abs(normW(realind,1:Ncan)),2),'b');
                hold on;
                h2=plot(usefreq,abs(normW(realind,1:Ncan)));

                title(sprintf('Sum %d modes, cord',Ncan))
            end; % else
        end; % if
    end;  % for orientmode

end; % if COMB3AXES

end; % for filenames



