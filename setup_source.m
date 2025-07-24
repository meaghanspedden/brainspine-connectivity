%% Setup for beamforming

sub='OP00219';
analysis='merged';

addpath('D:\for_meaghan\')
addpath('D:\msg_coreg')
addpath('D:\scannercast\tableofinfo') %for the stlread function


save_dir = fullfile('D:\MSST001', [sub '_' analysis]);
filename = fullfile('D:\MSST001', ...
    ['sub-' sub], ...
    'ses-001', ...
    'meg', ...
    ['p' analysis 'oe1000mspddfflo45hi45hfcstatic_001_array1.mat']);

%% get sensor positions

D=spm_eeg_load(filename);
grad=D.sensors('MEG');

%% Load and process the optical scan
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
% is in opm sensor space
fprintf('Select L shoulder R shoulder chin\n')
sim_fids_select = spm_mesh_select(mesh); %select left shoulder, right shoulder and chin
sim_fids = sim_fids_select';

S = [];
S.subject = mesh; 
S.fiducials = sim_fids;
S.plot = 'true';

T = cr_register_torso(S);

S = []; 
S.subject = mesh;
S.T = T;
S.sensors = grad;
S.spine_mode = 'default';
all_meshes = cr_check_registration(S);

%now we can save individual meshes in the space of opm sensor space
mesh_torso = all_meshes.mesh_torso;
mesh_wm = all_meshes.mesh_spine;
mesh_bone = all_meshes.mesh_vertebrae;
mesh_heart = all_meshes.mesh_heart;
mesh_lungs = all_meshes.mesh_lungs;

%create a source grid along the centerline of the spinal cord

%y_min = min(mesh_wm.vertices(:,2));
y_max = max(mesh_wm.vertices(:,2));
y_min=-32;


S = [];
S.spine = mesh_wm;
S.T = T;
S.resolution = 2.5; 
S.ylim = [y_min y_max];
S.unit = 'mm';
sources_center_line = cr_generate_spine_center(S);

figure; ft_plot_mesh(mesh_torso, 'facealpha',0.2, 'edgealpha', 0.2);hold on
plot3(sources_center_line.pos(:,1), sources_center_line.pos(:,2), sources_center_line.pos(:,3),'LineWidth',2)
axis image

save(fullfile(save_dir,'geoms'), 'T', 'sources_center_line', 'mesh_torso',...
    'mesh_wm', 'mesh_bone', 'mesh_heart', 'mesh_lungs')

%from here we now have a transform matrix (T), source grid and all relevent meshes in opm sensor space - this is everything we need for BEM and FEM forward modelling!