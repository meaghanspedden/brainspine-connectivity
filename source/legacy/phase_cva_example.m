close all;clear all;

%need to add SPM to path
spmpath='D:\spm';
addpath(spmpath)

N=1000; %% number trials
Nt=500; %% number time points per trial
Nf=Nt/20; %
Fs=1000; %% sample rate
T0=40; 
T1=20;
t=1:Nt; %% samples
x=zeros(N,Nt);y=x;
noise_range=[0:3];

figure;
for n1=1:length(noise_range) %% noisy sinnusoid amplitude
    UNITAMP=0; %% keeps constant magnitude within trials if set

    for tr=1:N

        a1=randn(1); %% random amplitudes
        a2=noise_range(n1)*randn(1); %% size of noisy sinusoild
        a3=randn(1); %% size of envelope modulation

        phase1=rand(1)*2*pi; %% also make
        phase2=rand(1)*2*pi;
        %% phase locked                           interferer                   
        y(tr,:)=a1*cos(2*pi*t/T0)+randn(size(t))+a2*sin(2*pi*t/T0)+...
            a3*sin(2*pi*t/T1+phase2); %% envelope covariation
        x(tr,:)=a1*sin(2*pi*t/T0)+randn(size(t))+...
            a3*sin(2*pi*t/T1+phase1); %% envelope covariation
        



        x(tr,:)=x(tr,:)-mean(x(tr,:));
        y(tr,:)=y(tr,:)-mean(y(tr,:));
        if UNITAMP,
            x(tr,:)=x(tr,:)./sqrt(dot(x(tr,:),x(tr,:)));
            y(tr,:)=y(tr,:)./sqrt(dot(y(tr,:),y(tr,:)));
        end

        fy=fft(y(tr,:));
        fx=fft(x(tr,:));

        fY(tr,:)=fy(2:Nf);
        fX(tr,:)=fx(2:Nf);
    end


    for f=1:Nf-1
        %% Now processing one freq band at a time

        %% remove DC level
        fY(:,f)=fY(:,f)-mean(fY(:,f));
        fX(:,f)=fX(:,f)-mean(fX(:,f));

        %% coherence
        Gxy=sum(fX(:,f).*conj(fY(:,f)));
        Gxx=sum(fX(:,f).*conj(fX(:,f)));
        Gyy=sum(fY(:,f).*conj(fY(:,f)));
        coh(f)=abs(Gxy).^2/(Gxx.*Gyy);
        %% regression coeff at each freq (similar to coherence)
        [B,BINT,R,RINT,STATS] = regress(fY(:,f),[fX(:,f) ones(size(fX(:,f)))]);
        Yrg=R;
        preg(f)=STATS(3);
        R2rg(f)=STATS(1);

        %% CVA
        %% split real and complex parts into columns of X and Y
        X=[real(fX(:,f)) imag(fX(:,f))];
        Y=[real(fY(:,f)) imag(fY(:,f))];
        %% get any linear mapping between X and Y
        CVA=spm_cva(Y,X);
        cvap(f)=CVA.p(1); % just look at first variate
        chicva(f,:)=CVA.chi;
        allCVAr(f,:)=CVA.r;
        v1=X*CVA.W(:,1); %% canonical variates
        v2=Y*CVA.V(:,1);
        R=corr(v1,v2); %% check correlation of variates
        r2x(f)=R.^2; %% this should be CVA.r.^2
        Nsig(f)=length(find(CVA.p<0.05)); %% number significant components
        r2_lin(f)=CVA.r(1).^2; %% linear variance explained


        cvashift=CVA.W*pinv(CVA.V); %% prediction coeffs
        Ypred=X*cvashift; %% prediction of Y based on X

        Yr=Y-Ypred; %% residuals after best linear prediction removed

        CVAr=spm_cva(Yr,X); %% double check- should be nothing significant here
        if CVAr.chi~=0
            error('CVA resids are correlated')
        end
        %% have removed linear phase-locked components,  now look for envelope
        Ye=abs([Yr(:,1)+Yr(:,2)*i]); % envelope
        Xe=abs([X(:,1)+X(:,2)*i]);
        Xe=Xe-mean(Xe);
        Ye=Ye-mean(Ye)

        [B,BINT,R,RINT,STATS] =regress(Ye,[Xe ones(size(Xe))]);
        r2env(f)=STATS(1) ;

    end % for f


    
    
    %h=plot(2:Nf,r2_lin,'x-',2:Nf,coh,'-o',2:Nf,r2env,':'); %% could plot envelope if needed
    subplot(2,2,n1)
    h=plot(2:Nf,r2_lin,'-',2:Nf,coh,':');
    legend('CVA','Coh')
    set(gca,'Fontsize',12)
    set(h,'LineWidth',3)
    set(h,'MarkerSize',12)
    title(sprintf('Noise=%3.2f',noise_range(n1)));
    xlabel('Frequency')
    ylabel('R^2')
end % for noise

