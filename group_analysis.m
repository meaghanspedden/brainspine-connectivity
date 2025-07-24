% group analysis

subs={'OP00212', 'OP00213', 'OP00214', 'OP00215', 'OP00219','OP00220', 'OP00221'};
which_ori='1';

t_vals = nan(length(subs),52);  % assuming each subject has a 1x40 vector of t-values

for i = 1:length(subs)
    sub = subs{i};
    save_dir = fullfile('D:\MSST001', [sub '_contrast']);
    savename = sprintf('spine_source_%s_%s.mat', sub, which_ori);
    
    % Load the data
    loaded = load(fullfile(save_dir,savename));  % adjust if variables are nested
    
    t_vals(i, 1:length(loaded.tvals)) = loaded.tvals;  % assuming size [1 x 40]

end

[n_participants, n_points] = size(t_vals);
[max_t_vals, peak_indices] = max(t_vals, [], 2);  % max and index per row

bin_width = 5;
edges = 1:bin_width:(n_points + 1);  % edges for bins of width 3

histogram(peak_indices, edges);
xlabel('Spinal Cord Source Point (binned)');
ylabel('Number of Participants');
title(sprintf('Peak t-value Locations (bin width = %d)', bin_width));

%%

peak_indices = sort(peak_indices);  % ensure sorted
cluster_radius = 3;  % hw many source points apart can still be considered a group
clusters = {};  % to store groups
used = false(size(peak_indices));  % logical index for processed peaks

for i = 1:length(peak_indices)
    if used(i)
        continue;
    end
    
    ref = peak_indices(i);
    
    % Group all peaks within cluster_radius of the reference
    close_idx = abs(peak_indices - ref) <= cluster_radius;
    
    % Store this group
    clusters{end+1} = peak_indices(close_idx);
    
    % Mark them as used
    used(close_idx) = true;
end



%% peak value density plot
y_coords = src.pos(peak_indices, 2);

% Compute KDE 
[f, xi] = ksdensity(y_coords, 'Bandwidth', 3);  % try 5cm, adjust as needed

% Plot
figure;
plot(xi, f, 'color',[0.24, 0.70, 0.44], 'LineWidth', 2); hold on;

% Overlay raw points at y-coordinate locations
scatter(y_coords, zeros(size(y_coords)), 100, 'v', 'filled', 'MarkerFaceColor', [0.49, 0.18, 0.56]);

xlabel('Spinal Cord Position (cm)');
ylabel('Density');
box off


%save(fullfile(save_dir,savename), 'pSig', 'tvals','src','vecs','thresholds','peakind')