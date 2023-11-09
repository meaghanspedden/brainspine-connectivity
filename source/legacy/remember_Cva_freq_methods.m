
fY=Ybrain_complex;
fX= emgdat_complex;
R2rg=zeros(1,numel(usefreq));
r2x=zeros(1,numel(usefreq));
r2_lin=zeros(1,numel(usefreq));
Nsig=zeros(1,numel(usefreq));
r2_env_CVA=zeros(1,numel(usefreq));
Nsigenv=zeros(1,numel(usefreq));
r2env=zeros(1,numel(usefreq));
resid_reg=zeros(size(fY));
resid_CVA=zeros(size(fY));
r2env_2=zeros(1,numel(usefreq));

for f=1:numel(usefreq)
    %% look for phase relationship first
    fY(:,f)=fY(:,f)-mean(fY(:,f));
    fX(:,f)=fX(:,f)-mean(fX(:,f));
    [B,BINT,Res,RINT,STATS] = regress(fY(:,f),[fX(:,f) ones(size(fX(:,f)))]);
    resid_reg(:,f)=Res;
    R2rg(f)=STATS(1);

    X=[real(fX(:,f)) imag(fX(:,f))];
    Y=[real(fY(:,f)) imag(fY(:,f))];
    X=X-mean(X);
    Y=Y-mean(Y);
    
    CVA=spm_cva(Y,X);

    v1=X*CVA.W(:,1); %% canonical variates
    v2=Y*CVA.V(:,1);

    R=corr(v1,v2); %% check correlation of variates
    r2x(f)=R.^2; %% this should be CVA.r.^2
    Nsig(f)=length(find(CVA.p<0.05)); %% number significant components
    r2_lin(f)=CVA.r(1).^2; %% linear variance explained


    cvashift=CVA.W/CVA.V; %% prediction coeffs
    Ypred=X*cvashift; %% prediction of Y based on X
    Yr=Y-Ypred; %% residuals after best linear prediction removed
    CVAr=spm_cva(Yr,X); %% double check- should be nothing significant here
    resid_CVA(:,f)=[Yr(:,1)+Yr(:,2)*i];

    %% now got rid of linear, look for envelope
    Ye=abs([Yr(:,1)+Yr(:,2)*i]); % envelope
    Xe=abs([X(:,1)+X(:,2)*i]);
    Xe=Xe-mean(Xe);
    Ye=Ye-mean(Ye);
    
    CVAenv=spm_cva(Ye,Xe); %% double check- should be nothing significant here
    Nsigenv(f)=length(find(CVAenv.p<0.05)); %% number significant components
    r2_env_CVA(f)=CVAenv.r(1).^2; %% linear variance explained


    cvashiftenv=CVAenv.W/CVAenv.V; %% prediction coeffs
    Ypredenv=Xe*cvashiftenv; %% prediction of Y based on X
    Yrenv=Ye-Ypredenv; %% residuals after best linear prediction removed


[B,BINT,R,RINT,STATS] =regress(Ye,[Xe ones(size(Xe))]);
r2env(f)=STATS(1) ;
%% do the same thing with theother residuals to compare
    Ye=abs(Res); % envelope
    Xe=abs([X(:,1)+X(:,2)*i]);
    Xe=Xe-mean(Xe);
    Ye=Ye-mean(Ye);
[B,BINT,~,RINT,STATS] =regress(Ye,[Xe ones(size(Xe))]);
r2env_2(f)=STATS(1) ;


end % for f
figure;
plot(usefreq,r2_lin,usefreq,r2env,usefreq,r2x,'x',usefreq,R2rg,'o')
legend('cva linear','envelope','can variates','regress')

figure
plot(usefreq,r2env_2,'ro-',usefreq,r2env,'b-')
legend({'complex reg residuals','CVA residuals'})
title('envelope correlations using residuals from each linear method')



