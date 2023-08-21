%error stats

subjectIDs={'116','122','123'};
meanerrorsR=[]; meanerrorsL=[];
stderrorsR=[]; stderrorsL=[];
for k=1:length(subjectIDs)
    subjectID=subjectIDs{k};
    savepath=['D:\MSST001\Coh_results00',subjectID];
    a=load(fullfile(savepath,sprintf('%s_error_all_r',subjectID))); %load trial wise errors RIGHT
    meanerrorsR=[meanerrorsR mean(a.allErrors)];
    stderrorsR=[stderrorsR std(a.allErrors)];
    b=load(fullfile(savepath,sprintf('%s_error_all_l',subjectID)));
    meanerrorsL=[meanerrorsL mean(b.allErrors)];
    stderrorsL=[stderrorsL std(b.allErrors)];

end
% Create a bar plot
x = 1:length(subjectIDs);
barWidth = 0.35;

figure;
hold on;

% Plot mean errors
bar(x - barWidth/2, meanerrorsR, barWidth, 'b', 'DisplayName', 'Right');
bar(x + barWidth/2, meanerrorsL, barWidth, 'r', 'DisplayName', 'Left');

% Error bars
errorbar(x - barWidth/2, meanerrorsR, stderrorsR, 'k.', 'LineWidth', 1.5);
errorbar(x + barWidth/2, meanerrorsL, stderrorsL, 'k.', 'LineWidth', 1.5);

hold off;

xlabel('Subject IDs');
ylabel('Errors (AU)');
title('Error for Right and Left Conditions');
legend({'Right','Left'},'Location', 'best');
xticks(x);
xticklabels(subjectIDs);
xlim([0.5, length(subjectIDs) + 0.5]);

grid on;
