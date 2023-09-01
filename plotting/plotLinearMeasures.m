function plotLinearMeasures(usefreq, Fstat_real, Fstat_imag, realSig, imagSig, coh,lab1,lab2)

figtitle=sprintf('%s %s linear connectivity',lab1,lab2);

cols=colormap(brewermap([],"Dark2"));
col3=cols(3,:); col4=cols(4,:); %for real imag and coh
col5=cols(8,:);

% Where to put significance stars
[maxValue, ~] = max(max([Fstat_real;Fstat_imag]));
yValue1 = maxValue + 0.1 * maxValue; %
yValue2 = yValue1 + 0.1 * maxValue;

figure;
subplot(2, 1, 1);
plot(usefreq, Fstat_real, 'color',col3,'LineWidth',3,'LineStyle',':'); hold on
plot(usefreq, Fstat_imag, 'color',col4, 'LineWidth', 3);
if ~isempty(realSig)
    plot(usefreq(realSig), ones(1,length(realSig))*yValue1,'*','MarkerSize',8,'color',col3)
end
if ~isempty(imagSig)
    plot(usefreq(imagSig), ones(1,length(imagSig))*yValue2,'*','MarkerSize',8,'color',col4)
end
ylabel('F Statistic');
xlabel('Frequency (Hz)');
legend('Real', 'Imaginary','Location','best');
legend boxoff
ax = gca;
ax.FontSize = 14;
ax.LineWidth=1.5; %change to the desired value     
box off
xlim([0 40])

% Plot Coherence
subplot(2, 1, 2);
if size(coh.cohspctrm,1) > 1
    coh.cohspctrm=squeeze(coh.cohspctrm(1,2,:));
end
plot(usefreq, coh.cohspctrm, 'LineWidth', 3,'color',col5);
ylabel('Coherence');
xlabel('Frequency (Hz)');

ax = gca;
ax.FontSize = 14;
ax.LineWidth=1.5; %change to the desired value     
box off
xlim([0 40])
sgtitle(figtitle, 'FontSize', 16, 'FontWeight', 'bold');
set(gcf, 'Position', [100, 100, 600, 800]);  % Adjust figure size and position


end
