
%% bootstrap cva and make barplot with error bars

clear all;close all;

subjectID='116';
figsavedir='D:\FiguresForPaper';

addpath(genpath('D:\brainspineconnectivity\plotting'))

whichanalysis={'brainemg','cordemg','braincord'};

cols=colormap(brewermap([],"Paired"));
col1=cols(1,:); col2=cols(2,:); 

RIGHTHAND=1;
LEFTHAND=~RIGHTHAND;

if isequal(RIGHTHAND,1)
    whichhand='r';
else
    whichhand='l';
end

if LEFTHAND
    cvaemg_l=['D:\MSST001\Coh_results00',subjectID,'\pplhandoe1_emg+abs.mat'];
    cvabrain_l=['D:\MSST001\Coh_results00',subjectID,'\pplhandoe1_brainopt+abs.mat'];
end
if RIGHTHAND
    %load cva data
    cvaemg_l=['D:\MSST001\Coh_results00',subjectID,'\pprhandoe1_emg+abs.mat'];
    cvabrain_l= ['D:\MSST001\Coh_results00',subjectID,'\pprhandoe1_brainopt+abs.mat'];
    %load regression data
    brainemg_reg=['D:\MSST001\Coh_results00',subjectID,'\regdat_brainemg_',subjectID,'_emg+abs.mat'];
    emgcord_reg=['D:\MSST001\Coh_results00',subjectID,'\regdat_spine_',subjectID,'_brainopt+abs.mat'];
    cordbrain_reg=['D:\MSST001\Coh_results00',subjectID,'\regdat_spine_',subjectID,'_emg+abs.mat'];
end

%load the data
a=load(cvaemg_l);
CVAemg=a.allCVA{1}; %cord emg at point of peak power
CVAemgbrain=a.CVAemgbrain_abs; % y is brain here x is emg
b=load(cvabrain_l); %spinal cord and brain, x is cord, y is brain, x0 is orthog brain sig
CVAbrain=b.allCVA{1};

emgbrainreg=load(brainemg_reg);
emgcordreg=load(emgcord_reg);
braincordreg=load(cordbrain_reg);

%% first do CVA

for i=1:length(whichanalysis)

    if strcmp(whichanalysis{i},'brainemg')
        X=CVAemgbrain.X;
        Y=CVAemgbrain.Y;
        X0=CVAemgbrain.X0;
    elseif strcmp(whichanalysis{i},'cordemg')
        X=CVAemg.X;
        Y=CVAemg.Y;
        X0=CVAemg.X0;

    elseif strcmp(whichanalysis{i},'braincord')
        X=CVAbrain.X;
        Y=CVAbrain.Y;
        X0=CVAbrain.X0;

    end

    bootvals=[];
    ntrials=size(X,1);

    for k=1:100

        newdatIdx = randsample(ntrials,round(ntrials/2)) ; %with replacement gives
        newdatX=[];
        newdatY=[];
        newdatX0=[];

        for j=1:round(ntrials/2) %reshape the data
            newdatX(j,:)=X(newdatIdx(j),:);
            newdatY(j,:)=Y(newdatIdx(j),:);
            newdatX0(j,:)=X0(newdatIdx(j),:);
        end

        % run the CVA
        CVA_boot=spm_cva(newdatY,newdatX, newdatX0);
        %calc var exp
        Ncan=max(find(CVA_boot.p<0.05));

        bootvals(k)=sum(CVA_boot.r(1:Ncan).^2);
    end


    var_nonlin(i)=mean(bootvals);
    std_nonlin(i)=std(bootvals);

end

%% now do linear
var_lin=[];
std_lin=[];
for i=1:length(whichanalysis)

    if strcmp(whichanalysis{i},'brainemg')
        Y=emgbrainreg.Ybrain2;
        X=emgbrainreg.Zemg;
        cind=emgbrainreg.cind;
        rpcind=emgbrainreg.rpcind;
    elseif strcmp(whichanalysis{i},'cordemg')
        Y=emgcordreg.Y;
        X=emgcordreg.X;
        cind=emgcordreg.cind;
        rpcind=emgcordreg.rpcind;

    elseif strcmp(whichanalysis{i},'braincord')
        Y=braincordreg.Y;
        X=braincordreg.X;
        cind=braincordreg.cind;
        rpcind=braincordreg.rpcind;

    end

    bootvals_lin=[];
    ntrials=size(X,1);

    for k=1:100

        newdatIdx = randsample(ntrials,round(ntrials/2)) ; %with replacement gives
        newdatX=[];
        newdatY=[];

        for j=1:round(ntrials/2) %reshape the data
            newdatX(j,:)=X(newdatIdx(j),:);
            newdatY(j,:)=Y(newdatIdx(j),:);
        end

        % run the regression
        rYbrain=[];
        for f=1:size(cind,2)*2

            [B,BINT,R,RINT,STATS] = regress(newdatY(:,f),[newdatX(:,rpcind(1,f)) newdatX(:,rpcind(2,f)) ones(size(newdatX,1),1)]);
            rYbrain(:,f)=R;


        end % for f

        %calc var exp here
        bootvals_lin(k)=1-var(rYbrain(:))/var(Y(:)); %residuals as prop of orig brain signal

    end

var_lin(i)=mean(bootvals_lin);
std_lin(i)=std(bootvals_lin);

end




%%


f=figure;
pos=[680   312   560   666];

data_nonlinear=var_nonlin *100;
std_nonlinear=std_nonlin * 100;

data_linear=var_lin * 100;
std_linear=std_lin * 100;

%Create subplots
subplot(2,1,1)
bar_data = bar(data_linear);
title('Linear')
ylabel('Variance explained (%)')
set(gca,'Xticklabel',{'Brain-EMG', 'Spinal cord-EMG', 'Brain-spinal cord'});
set(bar_data, 'Facecolor', col2)
ax = gca;
ax.FontSize = 14;
ax.LineWidth = 1.5;
box off

%Add error bars
hold on
num_bars = numel(data_linear);
x = bar_data.XEndPoints;
% get rid ofthird argument if you want error bars up and down
errorbar(x, data_linear, zeros(size(data_linear)), std_linear, 'k', 'linestyle', 'none', 'linewidth', 1.5);
hold off

subplot(2,1,2)
bar_data = bar(data_nonlinear);
box off
ylabel('Variance explained (%)')
title('Nonlinear')

set(gca,'Xticklabel',{'Brain-EMG', 'Spinal cord-EMG', 'Brain-spinal cord'});

set(bar_data, 'Facecolor', col1)
ax = gca;
ax.FontSize = 14;
ax.LineWidth = 1.5;
box off

% Add error bars
hold on
num_bars = numel(data_nonlinear);
x = bar_data.XEndPoints;
errorbar(x, data_nonlinear, zeros(size(data_nonlinear)),std_nonlinear, 'k', 'linestyle', 'none', 'linewidth', 1.5);
hold off

% Set the figure position
f.Position = pos;

savename=sprintf('VE_sub%s_%s.pdf',subjectID,whichhand);
exportgraphics(gcf, fullfile(figsavedir,savename), 'Resolution', 600)


%%

% for f=1:Ncanb,
%     blabstr(f,:)=sprintf('vb%2d',f);
%     blabstr(f+Ncanb,:)=sprintf('wb%2d',f);
% end;

% Ncane=max(find(CVAemg.p<0.05));
% vemg=CVAemg.v(:,1:Ncane);
% wemg=CVAemg.w(:,1:Ncane);
% for f=1:Ncane,
%     elabstr(f,:)=sprintf('ve%2d',f);
%     elabstr(f+Ncane,:)=sprintf('we%2d',f);
% end;


%Ncaneb=max(find(CVAemgbrain.p<0.05));
% vemgbrain=CVAemgbrain.v(:,1:Ncaneb);
% wemgbrain=CVAemgbrain.w(:,1:Ncaneb);

% for f=1:Ncaneb
%     eblabstr(f,:)=sprintf('veb%d',f);
%     eblabstr(f+Ncaneb,:)=sprintf('web%d',f);
% end


% cmatrix=[vemg wemg vbrain wbrain vemgbrain wemgbrain]; %canonical variates
% figure;
% c1=corr(cmatrix);
% imagesc(c1.^2);colorbar;
% set(gca,'Xtick',1:size(c1,1));
% set(gca,'Xticklabel',[elabstr;blabstr;eblabstr]);
% set(gca,'Ytick',1:size(c1,1));
% set(gca,'Yticklabel',[elabstr;blabstr;eblabstr]);
% title('W=cord, e=emg, b=brain')

% varemgcord_nonlin=sum(CVAemg.r(1:Ncane).^2); %% sum of r square for significant modes
% varbraincord_nonlin=sum(CVAbrain.r(1:Ncanb).^2);
% varbrainemg_nonlin=sum(CVAemgbrain.r(1:Ncaneb).^2);

%  varemgcord_nonlin=sum(CVAemg.r(1:end).^2); %% sum of r square for all modes
%  varbraincord_nonlin=sum(CVAbrain.r(1:end).^2);
%  varbrainemg_nonlin=sum(CVAemgbrain.r(1:end).^2);

%bootstrap variance
%stats = bootstrp(100,@(x)[mean(x) std(x)],y);






% Define your data and standard deviations
% data_linear = [varemgcordlin varbraincordlin varemgbrainlin] * 100;
% std_linear = [std_emgcordlin std_braincordlin std_brainemglin] * 100;
%
% data_nonlinear = [varemgcord_nonlin varbraincord_nonlin varbrainemg_nonlin] * 100;
% std_nonlinear = [std_emgcord_nonlin std_braincord_nonlin std_brainemg_nonlin] * 100;

