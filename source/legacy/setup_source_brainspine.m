%% Example script to co-register between OPM sensor positions, Optical Scan, and Simulation mesh models

% This pipeline is implemented for a subject dataset where we have the subject scan already in OPM sensor
% space, and hence we can skip registering the subject mesh to the opm sensors and go straight to pulling in
% the mesh models

%%

%this just needs to be run once, use same co reg for all
sub='OP00212';
analysis='merged';

addpath('D:\spm')
spm('defaults','EEG')
addpath('D:\for_meaghan\')
addpath('D:\msg_coreg')
addpath('D:\brainspineconnectivity') %for the stlread function


save_dir = fullfile('D:\MSST001', ['generic_' analysis]);
if~exist(save_dir, 'dir')
    mkdir(save_dir)
end
filename = fullfile('D:\MSST001', ...
    ['sub-' sub], ...
    'ses-001', ...
    'meg', ...
    ['p' analysis 'oe1000mspddfflo45hi45hfcstatic_001_array1.mat']);

%% get sensor positions

D=spm_eeg_load(filename);
grad=D.sensors('MEG');

%% Load and process the optical scan
cd('D:\new_leadfields_and_geom')
optical = ft_read_headshape('surface.stl', 'unit', 'mm');
p = struct();
p.vertices = optical.pos;
p.faces = optical.tri;
p2 = reducepatch(p, 0.01); % Reduce complexity of the mesh
subject.pos = p2.vertices;
subject.tri = p2.faces;

% Prepare mesh structure
mesh = struct();
mesh.vertices = subject.pos;
mesh.faces = subject.tri;
mesh.unit = 'mm'; % Define unit explicitly - make sure it matches the grad structure!

%%
% now we need to register the simulation meshes to the optical scan which
% is in opm sensor space - but we know what the transform matrix is so can
% just go ahead and apply this
% 
S = []; 
S.subject = mesh;
S.sensors = sensors(D, 'MEG'); 
S.brain = true;
S.spine_mode = 'cervical'; %full or cervical
S.torso_mode = 'anatomical'; %anatomical or canonical (if choosing canonical - you will need to select three fiducials on the subject mesh (left shoulder, right shoulder, chin))
all_meshes = cr_check_registration(S);

%now we can save individual meshes in the space of opm sensor space
mesh_torso = all_meshes.torso;
mesh_spine = all_meshes.spine;
mesh_bone = all_meshes.bone;
mesh_heart = all_meshes.heart;
mesh_lungs = all_meshes.lungs;
mesh_brain = all_meshes.brain;
mesh_iskull = all_meshes.iskull;
mesh_oskull = all_meshes.oskull;

%% create a source grid along the centerline of the spinal cord

y_min = min(mesh_spine.vertices(:,2));
y_max = max(mesh_spine.vertices(:,2));

S = [];
S.spine = mesh_spine;
S.T = all_meshes.transform;
S.resolution = 5; 
S.ylim = [y_min y_max];
S.unit = 'mm';
sources_spine = cr_generate_spine_center(S);

%% create source grid in the brain

bb_min = min(mesh_brain.vertices);
bb_max = max(mesh_brain.vertices);

[x,y,z] = ndgrid(...
    linspace(bb_min(1), bb_max(1), 10), ...
    linspace(bb_min(2), bb_max(2), 10), ...
    linspace(bb_min(3), bb_max(3), 10));

cand = [x(:) y(:) z(:)];

inside = false(size(cand,1),1);
for i = 1:size(cand,1)
    inside(i) = tt_is_inside(cand(i,:), mesh_brain.vertices, mesh_brain.faces);
end

sources_brain = [];
sources_brain.pos = cand(inside,:);
sources_brain.inside = true(size(sources_brain.pos,1),1);
sources_brain.unit = 'm';

%% visualise all sources
figure;
hold on;

ft_plot_mesh(mesh_torso, 'facecolor', 'purple', 'facealpha', 0.2, 'edgecolor', 'none');
ft_plot_mesh(mesh_brain, 'facecolor', 'b', 'facealpha', 0.2, 'edgecolor', 'none');
ft_plot_mesh(mesh_spine, 'facecolor', 'r', 'facealpha', 0.2, 'edgecolor', 'none');

plot3(sources_spine.pos(:,1), sources_spine.pos(:,2), sources_spine.pos(:,3), ...
    'k.', 'MarkerSize', 10);

plot3(sources_brain.pos(:,1), sources_brain.pos(:,2), sources_brain.pos(:,3), ...
    'k.', 'MarkerSize', 10);

axis equal
xlabel('X (m)')
ylabel('Y (m)')
zlabel('Z (m)')
title('Torso, Brain, Cord meshes with sources')
view(3)
camlight; lighting gouraud;

hold off;
%% from here we now have a transform matrix (T), source grid and all relevent meshes in opm sensor space - this is everything we need for BEM and FEM forward modelling!

mesh_heart.vertices = mesh_heart.vertices; 
mesh_heart.faces = mesh_heart.faces;

mesh_bone.vertices = mesh_bone.vertices;
mesh_bone.faces = mesh_bone.faces;

mesh_lungs.vertices = mesh_lungs.vertices; 
mesh_lungs.faces = mesh_lungs.faces;

mesh_torso.vertices = mesh_torso.vertices;
mesh_torso.faces = mesh_torso.faces;

mesh_wm.vertices = mesh_spine.vertices;
mesh_wm.faces = mesh_spine.faces;

coils_3axis = sensors(D,'meg');

mesh_brain.vertices = mesh_brain.vertices;
mesh_brain.faces = mesh_brain.faces;

mesh_iskull.vertices = mesh_iskull.vertices;
mesh_iskull.faces = mesh_iskull.faces;

mesh_oskull.vertices = mesh_oskull.vertices;
mesh_oskull.faces = mesh_oskull.faces;

sources_cent = sources_spine;
sources_brain = sources_brain;% Save to .mat