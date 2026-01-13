
meshinfl=load("C:\Users\mspedden\Documents\fieldtrip\template\anatomy\surface_pial_both.mat");
mesh=meshinfl.mesh;

meshspm=[];
meshspm.pos=brain_gii.vertices;
meshspm.tri=brain_gii.faces;

T=load("C:\Users\mspedden\Documents\brainspine_save\brain_transform_matrix.mat"); %from MNI to native
T=T.T;
% mri=ft_read_mri('single_subj_T1.nii');
% 
% mesh.pos = mesh.pos ./ 1.5;

mesh_n = ft_transform_geometry(T, mesh);

figure; ft_plot_mesh(mesh,  'facecolor',[0.8 0.8 0.8], 'facealpha',0.4,'edgecolor','none');
camlight
hold on
ft_plot_mesh(meshspm, 'facecolor','r','facealpha',0.6,'edgecolor','none')

%% 

atlas=ft_read_atlas('ROI_MNI_V4.nii');

[x,y,z] = ndgrid(1:atlas.dim(1), 1:atlas.dim(2), 1:atlas.dim(3));
voxels = [x(:) y(:) z(:) ones(numel(x),1)]';  % homogeneous coordinates
mni_coords = atlas.transform * voxels;  % 4×N, last row = 1
mni_coords = mni_coords(1:3,:)';   

native_coords_h = (T * [mni_coords, ones(size(mni_coords,1),1)]')';  
native_coords = native_coords_h(:,1:3);  % Nx3 coordinates in native space

native_labels = zeros(size(sourcemodel.pos,1),1);

for i = 1:size(coh_source.pos,1)
    % find nearest voxel in native_coords
    [~, idx] = min(sum((native_coords - coh_source.pos(i,:)).^2,2));
    linear_idx = sub2ind(atlas.dim, x(idx), y(idx), z(idx));
    native_labels(i) = atlas.tissue(linear_idx);
end




% nearest-neighbor atlas label for each vertex
mesh_labels = zeros(size(mesh_brain.pos,1),1);

for v = 1:size(mesh_brain.pos,1)
    % find closest voxel in atlas
    distances = sum((native_coords - mesh_brain.pos(v,:)).^2,2);
    [~, idx] = min(distances);
    mesh_labels(v) = atlas.tissue(idx);
end