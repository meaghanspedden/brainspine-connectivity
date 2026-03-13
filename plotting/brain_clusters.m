%% plot_M1_roi_p1_overlap.m
% Plots M1 ROI sphere, Participant 1 significant sources, and their overlap.
% Loads saved Step 4 results — no need to re-run the pipeline.

%% ---- USER CONFIG --------------------------------------------------------
save_dir       = 'C:\Users\mspedden\Documents\brainspine_save_bemv2';
geomfile       = 'C:\Users\mspedden\Documents\new_leadfields_and_geom\geometries_cervical_realistic.mat';
fieldtrip_path = 'C:\Users\mspedden\Documents\fieldtrip';
T_mat_path     = 'C:\Users\mspedden\Documents\brainspine_save\T.mat';
result_file    = 'groupRes_brain_DICS_bemv2_brainEMG_brainSmooth_8mm.mat';

p1_idx         = 1;
M1_mni_centre  = [-38 -26 56];
M1_radius_mm   = 20;
dist_tol_mm    = 5;     % sources within this distance count as overlap
saveFigs       = true;
%% -----------------------------------------------------------------------

addpath(fieldtrip_path); ft_defaults;

load(fullfile(save_dir, result_file), 'subjResults');
load(geomfile); mesh_brain.unit = 'mm';
load(T_mat_path);

% Build M1 ROI sphere in native space
[xg,yg,zg]  = ndgrid(-M1_radius_mm:M1_radius_mm, ...
                      -M1_radius_mm:M1_radius_mm, ...
                      -M1_radius_mm:M1_radius_mm);
sphere_mask  = (xg.^2 + yg.^2 + zg.^2) <= M1_radius_mm^2;
roi_coords   = [xg(sphere_mask), yg(sphere_mask), zg(sphere_mask)] + M1_mni_centre;
roi_hom      = [roi_coords, ones(size(roi_coords,1),1)];
roi_native   = (T * roi_hom')'; roi_native = roi_native(:,1:3);

% Participant 1 significant sources
sr1       = subjResults(p1_idx);
sig_mask  = sr1.coh_diff > sr1.thr95;
sig_pos   = sr1.pos(sig_mask, :);
sig_pos   = sig_pos(~any(isnan(sig_pos),2), :);

% Overlap: P1 sig sources within dist_tol of ROI
if ~isempty(sig_pos)
    D           = pdist2(roi_native, sig_pos);
    roi_overlap = min(D,[],2) <= dist_tol_mm;
    p1_overlap  = min(D,[],1)' <= dist_tol_mm;
    overlap_pos = sig_pos(p1_overlap, :);
    fprintf('P1 significant sources: %d\n', size(sig_pos,1));
    fprintf('P1 sources within M1 ROI: %d / %d (%.0f%%)\n', ...
        sum(p1_overlap), size(sig_pos,1), 100*mean(p1_overlap));
else
    overlap_pos = [];
    fprintf('No significant sources for Participant 1.\n');
end

% Plot
hfig = figure('Color','w');
ft_plot_mesh(mesh_brain,'facecolor',[0.85 0.85 0.85],'facealpha',0.25,'edgecolor','none');
hold on;

% M1 ROI sphere
h1 = plot3(roi_native(:,1), roi_native(:,2), roi_native(:,3), 'o', ...
    'MarkerSize',3, 'MarkerFaceColor',[0.6 0.6 0.6], 'MarkerEdgeColor','none');

% P1 significant sources (non-overlapping)
if ~isempty(sig_pos)
    non_overlap_pos = sig_pos(~p1_overlap,:);
    h2 = plot3(non_overlap_pos(:,1), non_overlap_pos(:,2), non_overlap_pos(:,3), 'o', ...
        'MarkerSize',6, 'MarkerFaceColor',[0.2 0.4 0.8], 'MarkerEdgeColor','none');
end

% Overlap
if ~isempty(overlap_pos)
    h3 = plot3(overlap_pos(:,1), overlap_pos(:,2), overlap_pos(:,3), 'o', ...
        'MarkerSize',8, 'MarkerFaceColor',[0.1 0.8 0.2], 'MarkerEdgeColor','k','LineWidth',0.5);
end

axis equal; camlight; lighting gouraud; view(176,-10);
title(sprintf('P1 (%s) — M1 ROI overlap  |  %d sig. sources, %d in ROI', ...
    sr1.subjID, size(sig_pos,1), size(overlap_pos,1)), ...
    'Interpreter','none','FontSize',12);
legend({'M1 ROI sphere','P1 sig. sources','Overlap'}, ...
    'Location','best','FontSize',11,'Box','off');

fig_dir = fullfile(save_dir,'figures');
if saveFigs && ~exist(fig_dir,'dir'), mkdir(fig_dir); end
if saveFigs
    savefig(hfig, fullfile(fig_dir,'step4_P1_M1_ROI_overlap.fig'));
end

if ~isempty(overlap_pos)
    centroid_native = mean(overlap_pos, 1);
    T_inv = inv(T);
    h = [centroid_native, 1]';
    centroid_mni = round((T_inv * h)');
    fprintf('Centroid of P1-M1 overlap (native): [%.1f %.1f %.1f] mm\n', ...
        centroid_native(1), centroid_native(2), centroid_native(3));
    fprintf('Centroid of P1-M1 overlap (MNI)   : [%d %d %d]\n', ...
        centroid_mni(1), centroid_mni(2), centroid_mni(3));
end