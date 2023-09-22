function plotNonlinearMeasures(usefreq, V,W,lab1,lab2)

cols=colormap(brewermap([],"Dark2"));
col1=cols(1,:); col2=cols(2,:);

figure;
plot(usefreq,V,'LineWidth',3,'color',col1)
% xlabel('Frequency (Hz)')
% ylabel('Canonical vectors (AU)')
box off
ax = gca;
ax.FontSize = 18;
ax.LineWidth=1.5; %change to the desired value     
hold on
yyaxis right
plot(usefreq,W,'LineWidth',3,'color', col2,'LineStyle','-.')
%ylabel(lab2)
legend({lab1,lab2},'Location','best')
legend boxoff
xlim([0 40])

