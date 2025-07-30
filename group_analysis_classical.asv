

%loads output from source_freq_group (pow diff), smooths, and does t test across subjects

subs = {'OP00212', 'OP00213', 'OP00214', 'OP00215', 'OP00219', ...
    'OP00220', 'OP00221', 'OP00224', 'OP00225', 'OP00226'};
which_ori = 'all';

load('D:\MSST001\generic_merged\geoms.mat')
powdiffs=zeros(length(subs,50)); %50 is n sourcepoints

for i = 1:length(subs)
    sub = subs{i};
    save_dir = fullfile('D:\MSST001', [sub '_contrast']);
    savename = sprintf('spine_source_group%s_%s.mat', sub, which_ori);

    loaded = load(fullfile(save_dir, savename));
    span_mm = 35;  % desired smoothing in mm
    dy = median(diff(geoms.src.pos(:,2)));  % spacing in mm
    span_pts = round(span_mm / dy);  % convert mm to number of points

    smoothed_powdiffs = smoothdata(powdiffs, 2, 'movmean', span_pts);

    powdiffs(i,:)=smoothed_powdiffs;

end


[h, p, ci, stats] = ttest(powdiffs,0, 'Tail', 'right'); %only looking for pos changes (decrease in power)

figure
plot(src.pos(:,2), stats.tvals,'k')
sig_idx = find(p < 0.05);
scatter(y_coords(sig_idx), tvals(sig_idx), 50, 'k', 'LineWidth', 1.5, ...
    'MarkerFaceColor', 'none', 'MarkerEdgeColor', 'r');  % hollow circles
title('significant power decrease across subjects')






