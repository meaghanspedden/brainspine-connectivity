
clear all
close all
clc
%loads output from source_freq_group (pow diff), smooths, and does t test across subjects
%now uncorrected for mult comp

subs = {'OP00212', 'OP00213', 'OP00215', 'OP00219', ...
    'OP00220', 'OP00221', 'OP00224', 'OP00225', 'OP00226'};
which_ori = 'all';

HFC=1;
negpk=1; %look for pos or neg peaks

load('D:\MSST001\generic_merged\geoms.mat')
powdiffs=[]; %will be total number of trials by n sourcepoints (50)

for i = 1:length(subs)
    sub = subs{i};
    save_dir = fullfile('D:\MSST001', [sub '_contrast']);
    
    if HFC
        %savename = sprintf('spine_source_%s_%s.mat', sub, which_ori);
            savename = sprintf('spine_source_group_%s_%s_fixed.mat', sub, which_ori);

    else
        savename = sprintf('spine_source_group_%s_%s_nohfc_fixed.mat', sub, which_ori);

    end

    loaded = load(fullfile(save_dir, savename));
    span_mm = 2.5;  % desired smoothing in mm
    dy = abs(median(diff(sources_center_line.pos(:,2))));  % spacing in mm
    span_pts = abs(round(span_mm / dy));  % convert mm to number of points

         smoothed_powdiffs = smoothdata(loaded.powdiff, 2, 'gaussian', span_pts);
   
         powdiffs=[smoothed_powdiffs; powdiffs];

end

%% permutation

% log_powrest_all = []; log_powstat_all = [];
% 
% for i = 1:length(subs)
%     sub = subs{i}; save_dir = fullfile('D:\MSST001', [sub '_contrast']);
%     
%     if HFC
%         savename = sprintf('spine_source_%s_%s.mat', sub, which_ori);
%     else
%         savename = sprintf('spine_source_%s_%s_nohfc.mat', sub, which_ori);
%     end
%    
%     loaded = load(fullfile(save_dir, savename));
% 
%     span_mm = 30;
%     dy = abs(median(diff(sources_center_line.pos(:,2))));
%     span_pts = abs(round(span_mm / dy));

    % Smooth log powers 
    % smoothed_rest = smoothdata(loaded.rest, 2, 'gaussian', span_pts); 
   % smoothed_stat = smoothdata(loaded.static, 2, 'gaussian', span_pts);

%     log_powrest_all = [log_powrest_all; loaded.rest]; 
%     log_powstat_all = [log_powstat_all; loaded.static];
% end

% combined_all = [log_powrest_all; log_powstat_all];  % [2N_trials xsourcepoints] 
% labels = [zeros(size(log_powrest_all, 1), 1);
%           ones(size(log_powstat_all, 1), 1)];
% 
% n_permutations = 5000; num_trials = size(combined_all, 1); 
% num_sources =size(combined_all, 2); 
% true_diff = mean(log_powrest_all, 1) - mean(log_powstat_all, 1);  % observed %mean across all subjects at each sourcepoint is neg
% 
% null_diffs = zeros(n_permutations, num_sources);
% 
% for i = 1:n_permutations
%     shuffled_idx = randperm(num_trials); 
%     shuffled_labels =labels(shuffled_idx);
% 
%     group_static = combined_all(shuffled_labels == 1, :); 
%     group_rest = combined_all(shuffled_labels == 0, :);
% 
%     % Make sure groups are same size 
%     min_n = min(size(group_static, 1),size(group_rest, 1)); 
%     group_static = group_static(1:min_n, :);
%     group_rest = group_rest(1:min_n, :);
% 
%     null_diffs(i, :) = mean(group_rest, 1) - mean(group_static, 1);
% end
% 
% % Compute one-sided p-values (e.g., static > rest) 
% p_perm =mean(null_diffs >= true_diff, 1);
%

%% from paramtric stuff

if negpk
    [h, p, ci, stats] = ttest(powdiffs,0, 'Tail', 'right'); %only looking for pos changes (decrease in power)
else
    [h, p, ci, stats] = ttest(powdiffs,0, 'Tail', 'left'); %only looking for neg changes (increase in power)

end
ycoords=sources_center_line.pos(:,2);


% 
% figure;
% boxplot(powdiffs, 'PlotStyle', 'compact', 'Colors', 'k', 'Symbol', '.');
% xlabel('Sourcepoint Index');
% ylabel('Power Difference Across Trials');
% title('Distribution of Power Differences per Sourcepoint');
% grid on;

addpath('D:\Violinplot-Matlab')
violinplot(powdiffs);  % trials x sources
xlabel('Sourcepoint Index');
ylabel('Power Difference Across Trials');
title('Power diff distribution concat across subject ');

nBins = 50;
data=powdiffs;
allDataMin = min(data(:));
allDataMax = max(data(:));
binEdges = linspace(allDataMin, allDataMax, nBins + 1);

histMatrix = zeros(nBins, size(data, 2));

for i = 1:size(data, 2)
    histCounts = histcounts(data(:, i), binEdges);
    histMatrix(:, i) = histCounts';
end

% Plot heatmap
figure;
imagesc(histMatrix);
colormap('hot');
colorbar;
xlabel('Sourcepoint');
ylabel('Histogram Bin');
title('Histogram Heatmap of Each Sourcepoint');
binCenters = (binEdges(1:end-1) + binEdges(2:end)) / 2;
yticks = round(linspace(1, nBins, 5));  % e.g. 5 tick labels
yticklabels = arrayfun(@(x) sprintf('%.2f', binCenters(x)), yticks, 'UniformOutput', false);
set(gca, 'YTick', yticks, 'YTickLabel', yticklabels);
[~, zeroIdx] = min(abs(binCenters - 0));
hold on;
yline(zeroIdx, 'c-', 'LineWidth', 1.5);



figure
plot(sources_center_line.pos(:,2), stats.tstat,'k'); hold on
sig_idx = find(p < 0.05);
ylabel('t-stat')
scatter(ycoords(sig_idx), stats.tstat(sig_idx), 50, 'k', 'LineWidth', 1.5, ...
    'MarkerFaceColor', 'none', 'MarkerEdgeColor', 'r');  % hollow circles
title('power change across subjects')
xlabel('y position')
if ~negpk
    set(gca, 'YDir', 'reverse');
end







