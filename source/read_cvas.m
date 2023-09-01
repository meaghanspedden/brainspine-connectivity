clear all;close all;

subjectID='116';
figsavedir='D:\FiguresForPaper';


cols=colormap(brewermap([],"Paired"));
col1=cols(1,:); col2=cols(2,:); 

RIGHTHAND=1;
LEFTHAND=~RIGHTHAND;

if LEFTHAND
    cvaemg_l=['D:\MSST001\Coh_results00',subjectID,'\pplhandoe1_emg+abs.mat'];
    cvabrain_l=['D:\MSST001\Coh_results00',subjectID,'\pplhandoe1_brainopt+abs.mat'];
end
if RIGHTHAND
    cvaemg_l=['D:\MSST001\Coh_results00',subjectID,'\pprhandoe1_emg+abs.mat'];
    cvabrain_l= ['D:\MSST001\Coh_results00',subjectID,'\pprhandoe1_brainopt+abs.mat'];
end

a=load(cvaemg_l);
usefreq=a.usefreq; %freqs analyzed
CVAemg=a.allCVA{1}; %cord emg at point of peak power
CVAemgbrain=a.CVAemgbrain_abs;
Fstatlin_emgbrainri=a.Fstatlin_emgbrainri;
Fstatlin_emgcord=a.Fstatlinri(:,1);
varemgcordlin=a.linvar_rmvd;
varemgbrainlin=a.linvaremgbrain_rmvd;

b=load(cvabrain_l); %spinal cord and brain
CVAbrain=b.allCVA{1};
varbraincordlin=b.linvar_rmvd;
Fstatlin_braincord=b.Fstatlinri(:,1);

Ncanb=max(find(CVAbrain.p<0.05));

vbrain=CVAbrain.v(:,1:Ncanb);
wbrain=CVAbrain.w(:,1:Ncanb);
for f=1:Ncanb,
    blabstr(f,:)=sprintf('vb%2d',f);
    blabstr(f+Ncanb,:)=sprintf('wb%2d',f);
end;

Ncane=max(find(CVAemg.p<0.05));
vemg=CVAemg.v(:,1:Ncane);
wemg=CVAemg.w(:,1:Ncane);
for f=1:Ncane,
    elabstr(f,:)=sprintf('ve%2d',f);
    elabstr(f+Ncane,:)=sprintf('we%2d',f);
end;


Ncaneb=max(find(CVAemgbrain.p<0.05));
vemgbrain=CVAemgbrain.v(:,1:Ncaneb);
wemgbrain=CVAemgbrain.w(:,1:Ncaneb);

for f=1:Ncaneb
    eblabstr(f,:)=sprintf('veb%d',f);
    eblabstr(f+Ncaneb,:)=sprintf('web%d',f);
end


% cmatrix=[vemg wemg vbrain wbrain vemgbrain wemgbrain]; %canonical variates
% figure;
% c1=corr(cmatrix);
% imagesc(c1.^2);colorbar;
% set(gca,'Xtick',1:size(c1,1));
% set(gca,'Xticklabel',[elabstr;blabstr;eblabstr]);
% set(gca,'Ytick',1:size(c1,1));
% set(gca,'Yticklabel',[elabstr;blabstr;eblabstr]);
% title('W=cord, e=emg, b=brain')

varemgcord_nonlin=sum(CVAemg.r(1:Ncane).^2); %% sum of r square for significant modes
varbraincord_nonlin=sum(CVAbrain.r(1:Ncanb).^2);
varbrainemg_nonlin=sum(CVAemgbrain.r(1:Ncaneb).^2);

%  varemgcord_nonlin=sum(CVAemg.r(1:end).^2); %% sum of r square for all modes
%  varbraincord_nonlin=sum(CVAbrain.r(1:end).^2);
%  varbrainemg_nonlin=sum(CVAemgbrain.r(1:end).^2);

f=figure;
pos=[680   312   560   666];

subplot(2,1,1)
h=bar([varemgcordlin varbraincordlin varemgbrainlin].*100);
title('Linear')
ylabel('Variance explained (%)')
set(gca,'Xticklabel',strvcat('EMG-spinal cord','Brain-spinal cord','Brain-EMG'));
ylim([0 20])
set(h,'Facecolor',col2)
ax = gca;
ax.FontSize = 14;
ax.LineWidth=1.5; %c
box off

subplot(2,1,2)
h=bar([varemgcord_nonlin varbraincord_nonlin varbrainemg_nonlin].*100);
ylim([0 80])
box off
ylabel('Variance explained (%)')
title('Nonlinear')
set(gca,'Xticklabel',strvcat('EMG-spinal cord','Brain-spinal cord','Brain-EMG'));
set(h,'Facecolor',col1)
ax = gca;
ax.FontSize = 14;
ax.LineWidth=1.5; %c
f.Position=pos;

if RIGHTHAND
    whichhand='r';
else
    whichhand='l';
end

savename=sprintf('VE_sub%s_%s.pdf',subjectID,whichhand);
exportgraphics(gcf, fullfile(figsavedir,savename), 'Resolution', 600)



