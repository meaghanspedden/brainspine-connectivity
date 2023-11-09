function new_CVA_FC1(Ydat,Xdat,refdat,usefreq,whichanalysis,figsavedir,subjectID,cnd)


fY=Ydat;
fX=Xdat;


r2_lin=zeros(1,numel(usefreq));
r2_env_CVA=zeros(1,numel(usefreq));


Nsig=zeros(1,numel(usefreq));
Nsigenv=zeros(1,numel(usefreq));

resid_CVA=zeros(size(fY));

p_env=zeros(1,numel(usefreq));
p_lin=zeros(1,numel(usefreq));

for f=1:numel(usefreq)
    %% look for phase relationship first
    fY(:,f)=fY(:,f)-mean(fY(:,f));
    fX(:,f)=fX(:,f)-mean(fX(:,f));
    refdat(:,:,f)=refdat(:,:,f)-mean(refdat(:,:,f)); %is this mean over correct dimension

    X=[real(fX(:,f)) imag(fX(:,f))];
    Y=[real(fY(:,f)) imag(fY(:,f))];
    X0=[real(refdat(:,:,f)) imag(refdat(:,:,f))];

    X=X-mean(X);
    Y=Y-mean(Y);
    X0=X0-mean(X0);

    CVA=spm_cva(Y,X,X0);

    Nsig(f)=length(find(CVA.p<0.05)); %% number significant components
    r2_lin(f)=CVA.r(1).^2; %% linear variance explained
    p_lin(f)=CVA.p(1); % p value

    cvashift=CVA.W/CVA.V; %% prediction coeffs
    Ypred=X*cvashift; %% prediction of Y based on X
    Yr=Y-Ypred; %% residuals after best linear prediction removed

    %% now got rid of linear, look for envelope
    Ye=abs([Yr(:,1)+Yr(:,2)*i]); % envelope
    Xe=abs([X(:,1)+X(:,2)*i]);
    X0e=abs(refdat(:,:,f));
    
    Xe=Xe-mean(Xe);
    Ye=Ye-mean(Ye);
    X0e=X0e-mean(X0e);

    CVAenv=spm_cva(Ye,Xe,X0e);
    
    Nsigenv(f)=length(find(CVAenv.p<0.05)); %% number significant components
    r2_env_CVA(f)=CVAenv.r(1).^2; %% linear variance explained
    p_env(f)=CVAenv.p(1);

    cvashiftenv=CVAenv.W/CVAenv.V; %% prediction coeffs
    Ypredenv=Xe*cvashiftenv; %% prediction of Y based on X
    Yrenv=Ye-Ypredenv; %% residuals after best linear prediction removed

    resid_CVA(:,f)=Yrenv;



end % for frequency loop

%% final CVA abs

Ycross=resid_CVA-mean(resid_CVA); %already abs; now mean centre
Xcross=abs(Xdat); %take the envelope
Xcross=Xcross-mean(Xcross);
CVAcrossfreq=spm_cva(Ycross,Xcross);

fprintf('CVA cross freq env p = %.3f chi=%.3f\n',CVAcrossfreq.p(1),CVAcrossfreq.chi(1))


%% run p vals through FDR


[~, crit_p, ~, adj_p]=fdr_bh(p_lin(5:end),0.05,'dep','no');
linSig=find(adj_p<0.05);

[~, crit_p, ~, adj_p]=fdr_bh(p_env(5:end),0.05,'dep','no');
envSig=find(adj_p<0.05);



%% plot phase amplitude association
cols=colormap(brewermap([],"Dark2"));
delete(gcf) %colormap command opens a figure, close it

col3=cols(3,:); col4=cols(4,:);
col5=cols(8,:);

% Where to put significance stars
[maxValue, ~] = max(r2_lin);
yValue1 = maxValue + 0.1 * maxValue; %

figure; plot(usefreq,r2_lin,'color',col4,'LineWidth',3)
xlabel('Frequency (Hz)');
if ~isempty(linSig)
    hold on
    plot(usefreq(linSig+4), ones(1,length(linSig))*yValue1,'*','MarkerSize',8,'color',col4)
end
ax = gca;
ax.FontSize = 18;
ax.LineWidth=.5; 
box off
xlim([0 40])
ylabel('R^2')
title('Phase and amplitude')


savename=sprintf('PA_sub%s_%s_%g.pdf',subjectID,whichanalysis,cnd);
exportgraphics(gcf, fullfile(figsavedir,savename), 'Resolution', 600)


%% plot envelope association
[maxValue, ~] = max(r2_env_CVA);
yValue = maxValue + 0.1 * maxValue; %
figure; plot(usefreq,r2_env_CVA,'color',col5,'LineWidth',3)
if ~isempty(envSig)
    hold on
    plot(usefreq(envSig+4), ones(1,length(envSig))*yValue,'*','MarkerSize',8,'color',col5)
end

xlabel('Frequency (Hz)');
ax = gca;
ax.FontSize = 18;
ax.LineWidth=.5; %change to the desired value
box off
xlim([0 40])
ylabel('R^2')
title('Envelope')

savename=sprintf('ENV_sub%s_%s_%g.pdf',subjectID,whichanalysis,cnd);
exportgraphics(gcf, fullfile(figsavedir,savename), 'Resolution', 600)

end %function





