clear all;close all;
RIGHTHAND=1;
LEFTHAND=~RIGHTHAND;

if LEFTHAND
    cvaemg_l=[];
    cvabrain_l=[];
end
if RIGHTHAND
    cvaemg_l='D:\MSST001\Coh_results00122\pprhandoe1_emg+abs.mat';
    cvabrain_l= 'D:\MSST001\Coh_results00122\pprhandoe1_brainopt+abs.mat';
end

a=load(cvaemg_l);
usefreq=a.usefreq;
CVAemg=a.allCVA{1};
CVAemgbrain=a.CVAemgbrain_abs;
Fstatlin_emgbrainri=a.Fstatlin_emgbrainri;
Fstatlin_emgcord=a.Fstatlinri(:,1);
varemgcordlin=a.linvar_rmvd;
varemgbrainlin=a.linvaremgbrain_rmvd;

b=load(cvabrain_l);
CVAbrain=b.allCVA{1};
varbraincordlin=b.linvar_rmvd;
Fstatlin_braincord=b.Fstatlinri(:,1);

Ncanb=max(find(CVAbrain.p<0.05));
normV_brain=(cov(CVAbrain.Y))*CVAbrain.V(:,1:Ncanb);
normW_brain=(cov(CVAbrain.X))*CVAbrain.W(:,1:Ncanb);
normW_brain=normW_brain./max(max(abs(normW_brain)));
normV_brain=normV_brain./max(max(abs(normV_brain)));
vbrain=CVAbrain.v(:,1:Ncanb);
wbrain=CVAbrain.w(:,1:Ncanb);
for f=1:Ncanb,
    blabstr(f,:)=sprintf('vb%2d',f);
    blabstr(f+Ncanb,:)=sprintf('wb%2d',f);
end;

Ncane=max(find(CVAemg.p<0.05));
normV_emg=(cov(CVAemg.Y))*CVAemg.V(:,1:Ncane);
normV_emg=normV_emg./max(max(abs(normV_emg)));
normW_emg=(cov(CVAemg.X))*CVAemg.W(:,1:Ncane);
normW_emg=normW_emg./max(max(abs(normW_emg)));
vemg=CVAemg.v(:,1:Ncane);
wemg=CVAemg.w(:,1:Ncane);
for f=1:Ncane,
    elabstr(f,:)=sprintf('ve%2d',f);
    elabstr(f+Ncane,:)=sprintf('we%2d',f);
end;

% figure;
% subplot(2,1,1);
% plot(1:length(normV_emg),normV_emg,1:length(normW_emg),normW_emg);
% subplot(2,1,2);
% plot(1:length(normV_brain),normV_brain,1:length(normW_brain),normW_brain);


% figure;
% plot(1:length(normV_brain),normV_brain(:,1),1:length(normV_emg),normV_emg(:,1))
% %plot(1:length(normV_brain),CVAbrain.V(:,1),1:length(normV_emg),CVAemg.V(:,1))
% legend('V brain(1)','V emg(1)')
% figure;
% plot(1:length(normW_brain),normW_brain(:,1),1:length(normW_emg),normW_emg(:,1))
% legend('W brain(1)','W emg(1)')

Ncaneb=max(find(CVAemgbrain.p<0.05));
normV_emgbrain=(cov(CVAemgbrain.Y))*CVAemgbrain.V(:,1:Ncaneb);
normV_emgbrain=normV_emgbrain./max(max(abs(normV_emgbrain)));
normW_emgbrain=(cov(CVAemgbrain.X))*CVAemgbrain.W(:,1:Ncaneb);
normW_emgbrain=normW_emgbrain./max(max(abs(normW_emgbrain)));
vemgbrain=CVAemgbrain.v(:,1:Ncaneb);
wemgbrain=CVAemgbrain.w(:,1:Ncaneb);

for f=1:Ncaneb,
    eblabstr(f,:)=sprintf('veb%d',f);
    eblabstr(f+Ncaneb,:)=sprintf('web%d',f);
end;


% plot(usefreq,Fstatlin_emgcord,usefreq,Fstatlin_braincord,usefreq,Fstatlin_emgbrain);
% legend('emg-cord','brain-cord','emg-brain')
% title(' linear')
% figure;
% plot(usefreq,normV_emgbrain,usefreq,normW_emgbrain);
% title('emg-brain, non-lin')


cmatrix=[vemg wemg vbrain wbrain vemgbrain wemgbrain];
figure;
c1=corr(cmatrix);
imagesc(c1.^2);colorbar;
set(gca,'Xtick',1:size(c1,1));
set(gca,'Xticklabel',[elabstr;blabstr;eblabstr]);
set(gca,'Ytick',1:size(c1,1));
set(gca,'Yticklabel',[elabstr;blabstr;eblabstr]);
title('W=cord, e=emg, b=brain')

figure;
% varemgcord_nonlin=sum(CVAemg.r(1:Ncane).^2); %% sum of r square for significant modes
% varbraincord_nonlin=sum(CVAbrain.r(1:Ncanb).^2);
% varbrainemg_nonlin=sum(CVAemgbrain.r(1:Ncaneb).^2);
 varemgcord_nonlin=sum(CVAemg.r(1:end).^2); %% sum of r square for significant modes
 varbraincord_nonlin=sum(CVAbrain.r(1:end).^2);
 varbrainemg_nonlin=sum(CVAemgbrain.r(1:end).^2);

subplot(2,1,1)
h=bar([varemgcordlin varbraincordlin varemgbrainlin].*100);
title('linear')
ylabel('Var explained')
set(gca,'Xticklabel',strvcat('emg-cord','brain-cord','brain-emg'));

subplot(2,1,2)
h=bar([varemgcord_nonlin varbraincord_nonlin varbrainemg_nonlin].*100);

ylabel('Var explained')
title('non-linear')
set(gca,'Xticklabel',strvcat('emg-cord','brain-cord','brain-emg'));
set(h,'Facecolor','r')



