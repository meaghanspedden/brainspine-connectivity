function plotLinearMeasures(usefreq,  Fstat, Sig,lab1,lab2,type)


 %put this in a function at some point!
    cols=colormap(brewermap([],"Dark2"));
    col3=cols(3,:); col4=cols(4,:);
    col5=cols(8,:);

    % Where to put significance stars
    [maxValue, ~] = max(Fstat);
    yValue1 = maxValue + 0.1 * maxValue; %

    figure;
    plot(usefreq, Fstat, 'color',col3,'LineWidth',3); hold on
    if ~isempty(Sig)
        plot(usefreq(Sig), ones(1,length(Sig))*yValue1,'*','MarkerSize',8,'color',col3)
    end
    if strcmp(type,'coh')
        ylabel('Coh')
    else
    ylabel('F Statistic');
    end
    xlabel('Frequency (Hz)');
    ax = gca;
    ax.FontSize = 18;
    ax.LineWidth=1.5; %change to the desired value
    box off
    xlim([0 40])
    title(sprintf('%s and %s %s',lab1,lab2,type))
