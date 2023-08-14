function plotNonlinearMeasures(usefreq, V,W,lab1,lab2)

cols=colormap(brewermap([],"Dark2"));
col1=cols(1,:); col2=cols(6,:);

figure;
plot(usefreq,V,'LineWidth',3,'color',col1)
xlabel('Frequency (Hz)')
ylabel(lab1)
box off
ax = gca;
ax.FontSize = 14;
hold on
yyaxis right
plot(usefreq,W,'LineWidth',3,'color', col2)
ylabel(lab2)
legend({lab1,lab2})
xlim([0 40])

