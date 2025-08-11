
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
powdiffs=zeros(length(subs),50); %50 is n sourcepoints

for i = 1:length(subs)
    sub = subs{i};
    save_dir = fullfile('D:\MSST001', [sub '_contrast']);
    if HFC
    savename = sprintf('spine_source_group_%s_%s.mat', sub, which_ori);
    else
    savename = sprintf('spine_source_group_%s_%s_nohfc.mat', sub, which_ori);

    end
    loaded = load(fullfile(save_dir, savename));
    span_mm = 20;  % desired smoothing in mm
    dy = abs(median(diff(sources_center_line.pos(:,2))));  % spacing in mm
    span_pts = abs(round(span_mm / dy));  % convert mm to number of points

    %smoothed_powdiffs = smoothdata(loaded.powdiff, 2, 'gaussian', span_pts);

    powdiffs(i,:)=loaded.powdiff;

end

if negpk
[h, p, ci, stats] = ttest(powdiffs,0, 'Tail', 'right'); %only looking for pos changes (decrease in power)
else
[h, p, ci, stats] = ttest(powdiffs,0, 'Tail', 'left'); %only looking for neg changes (increase in power)

end
ycoords=sources_center_line.pos(:,2);

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

%%
min_vals = zeros(size(powdiffs,1),1);
min_locs = zeros(size(powdiffs,1),1);

for subj = 1:size(powdiffs,1)
    [min_vals(subj), min_locs(subj)] = max(powdiffs(subj,:));
end

% Test if these negative peaks are significantly less than zero:
[h, p, ci, stats] = ttest(min_vals, 0, 'Tail', 'right');


[p, h, stats] = signrank(min_vals,0, 'tail', 'right');  % peaks = 1 value per subject
disp(p)










% %% permutation
% 
% data = powdiffs;  % 10 x 50
% 
% nSubjects = size(data,1);
% nSourcepoints = size(data,2);
% nPermutations = 1000;
% alpha = 0.05;
% 
% % Step 1: Compute observed t-statistics (one-tailed right)
% [~, ~, ~, stats] = ttest(data, 0, 'Tail', 'right');
% t_obs = stats.tstat;  % 1 x 50
% 
% % Step 2: Define cluster-forming threshold (critical t for one-tailed test)
% t_crit = tinv(1-alpha, nSubjects-1);
% 
% % Step 3: Find clusters above threshold
% % We'll treat sourcepoints as adjacent if they are neighbors (1D)
% above_thresh = t_obs > t_crit;
% 
% % Find contiguous clusters
% clusters = bwconncomp(above_thresh); % requires Image Processing Toolbox
% if clusters.NumObjects > 0
%     cluster_sums = zeros(1, clusters.NumObjects);
%     for i = 1:clusters.NumObjects
%         cluster_sums(i) = sum(t_obs(clusters.PixelIdxList{i}));
%     end
%     max_cluster_sum_obs = max(cluster_sums);
% 
% else
%     % No clusters found above threshold
%     max_cluster_sum_obs = -Inf;  % or -Inf if you prefer
% end
% 
% 
% % Step 4: Permutation test by sign flipping subjects
% max_cluster_sum_perm = zeros(1, nPermutations);
% for perm = 1:nPermutations
%     % Random sign flip for each subject
%     signs = randi([0 1], nSubjects, 1)*2 - 1;  % +1 or -1
%     permuted_data = data .* signs;
% 
%     % t-test on permuted data
%     [~, ~, ~, stats_perm] = ttest(permuted_data, 0, 'Tail', 'right');
%     t_perm = stats_perm.tstat;
% 
%     % Threshold and find clusters
%     above_thresh_perm = t_perm > t_crit;
%     clusters_perm = bwconncomp(above_thresh_perm);
% 
%     % Calculate max cluster sum for this permutation
%     if clusters_perm.NumObjects > 0
%         cluster_sums_perm = zeros(1, clusters_perm.NumObjects);
%         for j = 1:clusters_perm.NumObjects
%             cluster_sums_perm(j) = sum(t_perm(clusters_perm.PixelIdxList{j}));
%         end
%         max_cluster_sum_perm(perm) = max(cluster_sums_perm);
%     else
%         max_cluster_sum_perm(perm) = 0;
%     end
% end
% 
% % Step 5: Calculate cluster p-value
% p_value = mean(max_cluster_sum_perm >= max_cluster_sum_obs);
% 
% fprintf('Cluster-level p-value: %.4f\n', p_value);
% 
% 
% 
% 
% 
% 
