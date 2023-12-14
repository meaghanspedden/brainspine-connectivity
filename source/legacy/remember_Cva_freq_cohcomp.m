%% grb git test
clear all;
close all;
addpath D:\spm

N=1000;
Nt=512;
t=1:Nt;
x=zeros(N,Nt);y=x;
Nf=30;
for tr=1:N
    phase1=rand(1)*2*pi;
    phase2=rand(1)*2*pi;
    phase3=rand(1)*2*pi;
    phase4=rand(1)*2*pi;

    a1=randn(1);
    a2=randn(1);
    a3=randn(1);
    a4=randn(1);

%% linear relationship at t/10

    x(tr,:)=a1*sin(t/10+phase1)+randn(size(t)); %+a2*cos(t/10+phase3);
    y(tr,:)=a1*cos(t/10+phase1)+randn(size(t))+a2*sin(t/10); % a3*cos(t/10+phase2)
    x(tr,:)=x(tr,:)-mean(x(tr,:));
    y(tr,:)=y(tr,:)-mean(y(tr,:));
    
    fy=fft(y(tr,:));
    fx=fft(x(tr,:));

    fY(tr,:)=fy(2:Nf);
    fX(tr,:)=fx(2:Nf);
end;
figure
plot(mean(abs(fX)),'r');
ylabel('raw data')
xlabel('freq')

for f=1:Nf-1,
    %% look for phase relationship first
    fY(:,f)=fY(:,f)-mean(fY(:,f));
    fX(:,f)=fX(:,f)-mean(fX(:,f));
    [B,BINT,R,RINT,STATS] = regress(fY(:,f),[fX(:,f) ones(size(fX(:,f)))]);
    Yrg=R;
    R2rg(f)=STATS(1);
    
    Gxy=sum(fX(:,f).*conj(fY(:,f)));
    Gxx=sum(fX(:,f).*conj(fX(:,f)));
    Gyy=sum(fY(:,f).*conj(fY(:,f)));

    coh(f)=abs(Gxy).^2/(Gxx.*Gyy);


    X=[real(fX(:,f)) imag(fX(:,f)) rand(size(fX(:,f)))];
    Y=[real(fY(:,f)) imag(fY(:,f))];
 %   Y=[real(fY(:,f))];
    X=X-mean(X);
    Y=Y-mean(Y);
    
    CVA=spm_cva(Y,X);
    allCVAr(f,:)=CVA.r;
    v1=X*CVA.W(:,1); %% canonical variates
    v2=Y*CVA.V(:,1);
    R=corr(v1,v2); %% check correlation of variates
    r2x(f)=R.^2; %% this should be CVA.r.^2
    Nsig(f)=length(find(CVA.p<0.05)); %% number significant components
    r2_lin(f)=CVA.r(1).^2; %% linear variance explained


    cvashift=CVA.W*pinv(CVA.V); %% prediction coeffs
    Ypred=X*cvashift; %% prediction of Y based on X
%     R3=corr(Ypred,Y)
%     r3x(f)=R3(1,1);
    Yr=Y-Ypred; %% residuals after best linear prediction removed

    CVAr=spm_cva(Yr,X); %% double check- should be nothing significant here
    if CVAr.chi~=0, 
        error('CVA resids are correlated')
    end;
    %% now got rid of linear, look for envelope
    Ye=abs([Yr(:,1)+Yr(:,2)*i]); % envelope
    Xe=abs([X(:,1)+X(:,2)*i]);
    Xe=Xe-mean(Xe);
    Ye=Ye-mean(Ye)

    [B,BINT,R,RINT,STATS] =regress(Ye,[Xe ones(size(Xe))]);
    r2env(f)=STATS(1) ;
    
end; % for f
% figure;
% plot(2:Nf,r2_lin,2:Nf,r2env,2:Nf,r2x,'x',2:Nf,R2rg,'o')
% legend('cva linear','envelope','can variates','regress')
% 
% figure;
% plot(2:Nf,r2x,'-',2:Nf,R2rg,':x',2:Nf,coh,'d')
% legend('can variates','regress','coh')


figure;
% Plot the data with a specific line width
plot(2:Nf, r2x, 's-', 'LineWidth', 2);
hold on;  % This keeps the current plot and adds the next one to it
plot(2:Nf, coh, 'd-', 'LineWidth', 2);
ylabel('R^2')
legend('CVA','Coherence')
ax=gca;
ax.LineWidth=1.25;
legend box off


