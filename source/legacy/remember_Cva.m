
clear all;
close all;
addpath D:\spm

lookforstr='envelope';
%lookforstr='phase/amp';
N=1000;
Nt=512;
t=1:Nt;
x=zeros(N,Nt);y=x;
Nf=30;
for tr=1:N,
    phase1=rand(1)*2*pi;
    phase2=rand(1)*2*pi;
    rphase1=rand(1)*2*pi;
    rphase2=rand(1)*2*pi;

    a1=randn(1);
    a2=randn(1);
    a3=randn(1);
    %% linear interation at lowest freq (t/10)
    %% non-linear coupling between lowest and medium freq (t/5)
    %% just envelope coupling at hightest freq (t/3)
    x(tr,:)=a1*sin(t/10+phase1)+a3*sin(t/3+rphase1);%+randn(1,Nt);
    y(tr,:)=a1*cos(t/10+phase1)+a1*cos(t/5+phase2)+a3*sin(t/3+rphase2);%+randn(1,Nt);
    %y(tr,:)=a3*sin(t/3+rphase2);%+randn(1,Nt);
    fy=fft(y(tr,:));
    fx=fft(x(tr,:));
    
    fY(tr,:)=fy(2:Nf);
    fX(tr,:)=fx(2:Nf);
end;

 figure;
subplot(2,1,1);
plot(t,x(1,:))
subplot(2,1,2);
plot(t,y(1,:))

 
figure
plot(mean(abs(fX)),'r');
hold on
plot(mean(abs(fY)),'b');

ylabel('raw data')
xlabel('freq')

switch lookforstr,
    case 'phase/amp'
        X=[real(fX) imag(fX)]; % phase/amp
        Y=[real(fY) imag(fY)];
    case 'envelope',
        X=abs(fX); % envelope
        Y=abs(fY);
    case 'phase'
        X=angle(fX);
        Y=angle(fY);
end;
X=X-mean(X);
Y=Y-mean(Y);

CVA=spm_cva(Y,X);

Nsig=length(find(CVA.p<0.05)); %% number significant components

for f=1:Nsig,
    normV=cov(CVA.Y)*CVA.V;
    normW=cov(CVA.X)*CVA.W;
    figure;
    try,
        subplot(2,1,1);
        figure
        plot(1:Nf-1,normW(1:Nf-1,f),'g',1:Nf-1,normW(Nf:2*(Nf-1),f),'r');
        legend('real','imag')
        ylabel('W (design)')
    hold on;
        plot(1:Nf-1,normV(1:Nf-1,f),'g',1:Nf-1,normV(Nf:2*(Nf-1),f),'r');
        ylabel('V (data)')
        xlabel('freq')
    catch
        subplot(2,1,1);
        figure
       plot(1:Nf-1,normV(1:Nf-1,f),'b-o','MarkerFaceColor','b');
        hold on
        plot(1:Nf-1,normW(1:Nf-1,f),'r-o','MarkerFaceColor','r');
        ylabel('Canonical vectors (AU)')
        hold on
        xlabel('Frequency (Hz)')
        ax = gca;
        ax.FontSize = 18;
    end;

    title(sprintf('mode %d, p<%3.2f',f,CVA.p(f)))
    %% CVA.V(:,1) etc are the canonical vectors for the Y
    % CVA.v is the canonical variate (=Y*CVA.V)

    %% CVA.W(:,1) etc are the canonical vectors for the X
    % CVA.w is the canonical variate (=X*CVA.W)
    %% for some reason CVA.W is not set of eigenvectors (but CVA.V is).. need to check


end;

