% calc RMS error for each subject

clear all; close all; clc

sub='122';
hand='l';
savedir='D:\MSST001\Coh_results00122';
EMGdir=['D:\OP00',sub,'_experiment\EMGfiles\'];

if strcmp(sub,'116')
    targets=[67 52];
    rightfiles={'07', '09','11','13'}; %for 116
    leftfiles={'23', '25', '27'};
    runsleft={'001','002','003'};
    runsright={'001','002','003','004'};
    subsave='116';


elseif strcmp(sub,'122')
    targets=[58 62];
    rightfiles={'12','14','16','18'};
    leftfiles={'20','22','24','26'};
    runsright={'002','003','004','005'};
    runsleft={'001','002','003','004'};
    subsave='122';

elseif strcmp(sub,'123')
    targets=[63 62];
    rightfiles={'08','10','12','14'};
    leftfiles={'16','18','20','22'};
    runsright={'001','003','004','005'};
    runsleft={'001','002','003','004'};
    subsave='124';



else
    error('no subject with that ID')

end

if strcmp(hand,'r')
    emgfiles=rightfiles;
    target=targets(1);
    runs=runsright;
else
    emgfiles=leftfiles;
    target=targets(2);
    runs=runsleft;
end

allErrors=[];
for k=1:length(emgfiles)

    emgfile=[EMGdir,sub,'_',emgfiles{k}];
    load(emgfile)

    if strcmp(hand,'r')
        dat=Rect_smooth_EMG_r;
        retainedtrialsfile=fullfile(savedir,sprintf('retainedtrialsOP00%srun%s_r',subsave,runs{k}));


    elseif strcmp(hand,'l')
        dat=Rect_smooth_EMG_l;
        retainedtrialsfile=fullfile(savedir,sprintf('retainedtrialsOP00%srun%s_l',subsave,runs{k}));

    end

    load(retainedtrialsfile)


    %% synchronize

    trigtime=0.05;
    sf=1/dat.interval;
    trigsamp=round(trigtime*sf);
    synchedEMG=dat.values(trigsamp:end);
    newtime=dat.times(trigsamp:end);
    newtime=newtime-newtime(1);

    %% plot

    targetline=ones(1,length(synchedEMG))*target;

    figure; plot(newtime, synchedEMG); hold on
    plot(newtime,targetline,'r--','LineWidth',1.5)


    %calc RMSE for each second of data
    start=1;
    targetlineEpoch=targetline(1:sf);
    nepochs=floor(length(synchedEMG)/sf);
    disp(nepochs)
    for j=1:nepochs %fix this
        thisdat=synchedEMG(start:start+sf-1);
        epochErrors(j)=sqrt(mean((thisdat-targetlineEpoch').^2));
        start=start+sf;
    end

    %delete same trials deleted from OPM analysis
    epochErrors=epochErrors(retain);
    allErrors=[allErrors epochErrors];
end



savename=sprintf('%s_error_all_%s',sub,hand);
save(fullfile(savedir,savename),'allErrors')



