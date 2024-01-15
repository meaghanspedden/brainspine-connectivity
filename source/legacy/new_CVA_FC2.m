function [chi_lin, chi_env, chi_cross, df_lin, df_env, df_cross]= new_CVA_FC2(Ydat,Xdat,refdat,Prec,usefreq,freqtestidx,whichanalysis,figsavedir,subjectID,cnd)


fY=Ydat;
fX=Xdat;

r2_lin=zeros(1,numel(usefreq)); % save r2 values
r2_env_CVA=zeros(1,numel(usefreq));


Nsig=zeros(1,numel(usefreq)); %number sig components
Nsigenv=zeros(1,numel(usefreq));

resid_CVA=zeros(size(fY)); %save residuals for cross frequency analysis

p_env=zeros(1,numel(usefreq)); %save p values
p_lin=zeros(1,numel(usefreq));

chi_env=zeros(1,numel(usefreq));
chi_lin=zeros(1,numel(usefreq));

df_env=zeros(1,numel(usefreq));
df_lin=zeros(1,numel(usefreq));

for f=1:numel(usefreq)

    %% look for phase and amplitude relationship first


    %separate into real and imaginary parts
    X=[real(fX(:,f)) imag(fX(:,f))];
    Y=[real(fY(:,f)) imag(fY(:,f))];
    X0=[real(refdat(:,:,f)) imag(refdat(:,:,f))];

    if ~isempty(Prec)
        X=[X Prec];
    end

    % mean centre
    X=X-mean(X);
    Y=Y-mean(Y);
    X0=X0-mean(X0);

    CVA=spm_cva(Y,X,X0);

    Nsig(f)=length(find(CVA.p<0.05)); %% number significant components
    r2_lin(f)=CVA.r(1).^2; %% linear variance explained
    p_lin(f)=CVA.p(1); % p value
    chi_lin(f)=CVA.chi(1);
    df_lin(f)=CVA.df(1);

    cvashift=CVA.W/CVA.V; %% prediction coeffs
    Ypred=X*cvashift; %% prediction of Y based on X
    Yr=Y-Ypred; %% residuals after best linear prediction removed

    %test that relationship is gone
    CVAresid=spm_cva(Yr,X,X0);

    if CVAresid.p < 0.05
        error('there is still a relationship p should not be sig!')
    end


    %% now look for envelope
    %here we use residuals from above analysis

    Ye=abs([Yr(:,1)+Yr(:,2)*i]); % put back into complex and take abs
    Xe=abs([X(:,1)+X(:,2)*i]);
    X0e=abs(refdat(:,:,f));

    if ~isempty(Prec)
        Xe=[Xe Prec];
    end

    Xe=Xe-mean(Xe); %mean centre
    Ye=Ye-mean(Ye);
    X0e=X0e-mean(X0e);

    CVAenv=spm_cva(Ye,Xe,X0e);

    Nsigenv(f)=length(find(CVAenv.p<0.05)); %% number significant components
    r2_env_CVA(f)=CVAenv.r(1).^2; %% linear variance explained
    p_env(f)=CVAenv.p(1);
    chi_env(f)=CVAenv.chi(1);
    df_env(f)=CVAenv.df(1);

    cvashiftenv=CVAenv.W/CVAenv.V; %% prediction coeffs
    Ypredenv=Xe*cvashiftenv; %% prediction of Y based on X
    Yrenv=Ye-Ypredenv; %% residuals after best linear prediction removed

    CVAenvresid=spm_cva(Yrenv,Xe,X0e);
    if CVAenvresid.p < 0.05
        error('there is still a relationship p should not be sig!')
    end

    resid_CVA(:,f)=Yrenv; %save residuals for cross-freq analysis



end % for frequency loop

%% final CVA - cross frequency amplitude envelope

Ycross=resid_CVA-mean(resid_CVA); %already abs; now mean centre

Xcross=abs(Xdat); %take the envelope

if~isempty(Prec)
    Xcross=[Xcross Prec];
end

Xcross=Xcross-mean(Xcross);


% select freqs we want to test
Ycross=Ycross(:,freqtestidx:end);
Xcross=Xcross(:,freqtestidx:end);

CVAcrossfreq=spm_cva(Ycross,Xcross);
chi_cross=CVAcrossfreq.chi(1);
df_cross=CVAcrossfreq.df(1);

fprintf('CVA cross freq env p = %.5f chi=%.3f\n',CVAcrossfreq.p(1),CVAcrossfreq.chi(1))


%% now run p vals from per frequency analyses through FDR

if isempty(Prec)

    [~, ~, ~, adj_p]=fdr_bh(p_lin(freqtestidx:end),0.05,'dep','no');
    linSig=find(adj_p<0.05);

    [~, ~, ~, adj_p]=fdr_bh(p_env(freqtestidx:end),0.05,'dep','no');
    envSig=find(adj_p<0.05);


    %% plot phase coupling

    cols=colormap(brewermap([],"Dark2"));
    delete(gcf) %colormap command opens a figure, close it

    col1=cols(1,:); col2=cols(2,:);
    col3=cols(3,:); col4=cols(4,:);
    col5=cols(8,:);

    freqshift=freqtestidx-usefreq(1); %index shift because testing from 5 Hz

    % Where to put significance stars
    % [maxValue, ~] = max(r2_lin);
    % yValue1 = maxValue + 0.1 * maxValue;

    figure; pp=plot(usefreq,r2_lin,'color',col4,'LineWidth',3);
    %xlabel('Frequency (Hz)');


    if ~isempty(linSig)

        linSig=linSig+freqshift;

        discontinuity_indices = find(diff(linSig) > 1);

        % If there are no discontinuities, there's only one part
        if isempty(discontinuity_indices)
            Sigparts = {linSig};
        else
            % Split the vector at the discontinuity indices
            Sigparts = mat2cell(linSig, 1, diff([0, discontinuity_indices, length(linSig)]));
        end

        hold on
        %plot(usefreq(linSig), ones(1,length(linSig))*yValue1,'*','MarkerSize',8,'color',col4)

        for k=1:length(Sigparts)
            fill([usefreq(Sigparts{k}), fliplr(usefreq(Sigparts{k}))], ...
                [zeros(1, numel(Sigparts{k})), fliplr(r2_lin(Sigparts{k}))], col4, 'FaceAlpha', 0.3,'EdgeColor',col4,'LineWidth',2);

        end


    end

    ax = gca;
    ax.FontSize = 18;
    ax.LineWidth=1.25;
    box off
    xlim([0 40])
    %ylabel('R^2')
    %title('Phase and amplitude')
    ax = ancestor(pp, 'axes');
    ax.YAxis.Exponent = 0;
    %ytickformat('%.0f')

    savename=sprintf('PA_sub%s_%s_%g.pdf',subjectID,whichanalysis,cnd);
    exportgraphics(gcf, fullfile(figsavedir,savename), 'Resolution', 600)


    %% plot envelope association

    [maxValue, ~] = max(r2_env_CVA);
    yValue = maxValue + 0.1 * maxValue; %where to put significance stars

    figure; plot(usefreq,r2_env_CVA,'color',col5,'LineWidth',3)

    if ~isempty(envSig)

        envSig=envSig+freqshift;

        discontinuity_indices = find(diff(envSig) > 1);

        % If there are no discontinuities, there's only one part
        if isempty(discontinuity_indices)
            Sigparts = {envSig};
        else
            % Split the vector at the discontinuity indices
            Sigparts = mat2cell(envSig, 1, diff([0, discontinuity_indices, length(envSig)]));
        end

        hold on
        %plot(usefreq(linSig), ones(1,length(linSig))*yValue1,'*','MarkerSize',8,'color',col4)

        for k=1:length(Sigparts)
            fill([usefreq(Sigparts{k}), fliplr(usefreq(Sigparts{k}))], ...
                [zeros(1, numel(Sigparts{k})), fliplr(r2_env_CVA(Sigparts{k}))], col5, 'FaceAlpha', 0.1,'EdgeColor',col5,'LineWidth',2);

        end


    end

    %xlabel('Frequency (Hz)');
    ax = gca;
    ax.FontSize = 18;
    ax.LineWidth=1.25; %change to the desired value
    box off
    xlim([0 40])
    % ylabel('R^2')
    % title('Envelope')


    savename=sprintf('ENV_sub%s_%s_%g.pdf',subjectID,whichanalysis,cnd);
    exportgraphics(gcf, fullfile(figsavedir,savename), 'Resolution', 600)

    %% CVA cross freq plots for supplementary material

    %normalization according to Haufe

    %CVAcrossfreq
    Nsig=length(find(CVAcrossfreq.p<0.05)); %% number significant components
    if Nsig >= 1

        normV=cov(CVAcrossfreq.Y)*CVAcrossfreq.V(:,1:Nsig)*inv(cov(CVAcrossfreq.v(:,1:Nsig)));
        normW=cov(CVAcrossfreq.X)*CVAcrossfreq.W(:,1:Nsig)*inv(cov(CVAcrossfreq.w(:,1:Nsig)));


        if strcmp(whichanalysis,'emgbrain')

            lab1='Brain';
            lab2='EMG';

        elseif strcmp(whichanalysis,'cordemg')

            lab1='EMG';
            lab2='Spinal cord';

        else

            lab1='Brain';
            lab2='Spinal cord';

        end

        figure;
        plot(usefreq(freqtestidx:end),abs(normV(:,1)),'LineWidth',3,'color',col1) %Y
        box off
        ax = gca;
        ax.FontSize = 18;
        ax.LineWidth=1.25; %change to the desired value
        hold on
        yyaxis right
        plot(usefreq(freqtestidx:end),abs(normW(:,1)),'LineWidth',3,'color', col2,'LineStyle','-.') %X
        legend({lab1,lab2},'Location','best')
        legend boxoff
        xlim([0 40])

        savename=sprintf('CANVECS_sub%s_%s_%g.pdf',subjectID,whichanalysis,cnd);
        exportgraphics(gcf, fullfile(figsavedir,savename), 'Resolution', 600)


    else
        fprintf('CVA not significant\n')

    end

end %if prec

end %function





