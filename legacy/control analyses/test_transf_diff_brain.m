
meshinfl=load("C:\Users\mspedden\Documents\fieldtrip\template\anatomy\surface_inflated_both.mat");
%mesh2=load("C:\Users\mspedden\Documents\fieldtrip\template\anatomy\surface_pial_both.mat");

%mesh=meshinfl.mesh;
mesh=meshinfl.mesh;

mesh_n = ft_transform_geometry(T, mesh);

figure; ft_plot_mesh(mesh_n,  'facecolor',[0.8 0.8 0.8], 'facealpha',0.4,'edgecolor','none');
camlight
hold on
ft_plot_mesh(mesh_brain, 'facecolor','r','facealpha',0.6,'edgecolor','none')

%% 
% Center of left M1 hand knob
mni_center = [-38 -26 56];

% Radius in mm
radius = 20;

% Generate sphere in MNI space (for plotting or voxel selection)
[xg, yg, zg] = ndgrid(-radius:radius, -radius:radius, -radius:radius);
sphere_mask = (xg.^2 + yg.^2 + zg.^2) <= radius^2;

% Coordinates of voxels inside sphere
roi_coords = [xg(sphere_mask), yg(sphere_mask), zg(sphere_mask)] + mni_center;

% Convert to homogeneous coordinates
roi_hom = [roi_coords, ones(size(roi_coords,1),1)];

load('C:\Users\mspedden\Documents\brainspine_save\T.mat')

% Apply transformation
roi_native_hom = (T * roi_hom')';

% Extract native X/Y/Z
roi_native = roi_native_hom(:,1:3);

figure; ft_plot_mesh(mesh_brain,  'facecolor',[0.8 0.8 0.8], 'facealpha',0.7,'edgecolor','none');
hold on

% Plot the M1 ROI voxels in native space
plot3(roi_native(:,1), roi_native(:,2), roi_native(:,3), ...
      'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r', 'LineWidth', 1.5);

view([180, -20]);
camlight('headlight');
lighting gouraud;
axis equal;
xlabel('X'); ylabel('Y'); zlabel('Z');
title('M1 ROI in native space');















% atlas=ft_read_atlas('ROI_MNI_V4.nii');
% 
% [x,y,z] = ndgrid(1:atlas.dim(1), 1:atlas.dim(2), 1:atlas.dim(3));
% voxels = [x(:) y(:) z(:) ones(numel(x),1)]';  % homogeneous coordinates
% mni_coords = atlas.transform * voxels;  % 4×N, last row = 1
% mni_coords = mni_coords(1:3,:)';   
% 
% native_coords_h = (T * [mni_coords, ones(size(mni_coords,1),1)]')';  
% native_coords = native_coords_h(:,1:3);  % Nx3 coordinates in native space
% 
% native_labels = zeros(size(sourcemodel.pos,1),1);
% 
% for i = 1:size(coh_source.pos,1)
%     % find nearest voxel in native_coords
%     [~, idx] = min(sum((native_coords - coh_source.pos(i,:)).^2,2));
%     linear_idx = sub2ind(atlas.dim, x(idx), y(idx), z(idx));
%     native_labels(i) = atlas.tissue(linear_idx);
% end
% 
% 
% 
% 
% % nearest-neighbor atlas label for each vertex
% mesh_labels = zeros(size(mesh_brain.pos,1),1);
% 
% for v = 1:size(mesh_brain.pos,1)
%     % find closest voxel in atlas
%     distances = sum((native_coords - mesh_brain.pos(v,:)).^2,2);
%     [~, idx] = min(distances);
%     mesh_labels(v) = atlas.tissue(idx);
% end