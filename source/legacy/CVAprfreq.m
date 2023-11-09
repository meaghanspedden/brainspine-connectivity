
Nsig=zeros(1,numel(usefreq));
r2x=zeros(1,numel(usefreq));
r2_lin=zeros(1,numel(usefreq));

for f=1:numel(usefreq)
    %% look for phase relationship first
    Xi=[real(X(:,f)) imag(X(:,f))];
    Yi=[real(Y(:,f)) imag(Y(:,f))];
    Xi=Xi-mean(Xi);
    Yi=Yi-mean(Yi);
    CVA=spm_cva(Yi,Xi);
    v1=Xi*CVA.W(:,1); %% canonical variates
    v2=Yi*CVA.V(:,1);
    R=corr(v1,v2); %% check correlation of variates
    r2x(f)=R.^2; %% this should be CVA.r.^2
    Nsig(f)=length(find(CVA.p<0.05)); %% number significant components
    r2_lin(f)=CVA.r(1).^2; %% linear variance explained
% 
% 
%     cvashift=CVA.W/CVA.V; %% prediction coeffs
%     Ypred=Xi*cvashift; %% prediction of Y based on X
%     Yr=Y-Ypred; %% residuals after best linear prediction removed
%     CVAr=spm_cva(Yr,X); %% double check- should be nothing significant here

    %% now got rid of linear, look for envelope
%     Ye=abs([Yr(:,1)+Yr(:,2)*i]); % envelope
%     Xe=abs([X(:,1)+X(:,2)*i]);
%     Xe=Xe-mean(Xe);
%     Ye=Ye-mean(Ye)
% 
%     [B,BINT,R,RINT,STATS] =regress(Ye,[Xe ones(size(Xe))]);
%     r2env(f)=STATS(1) ;
end; % for f