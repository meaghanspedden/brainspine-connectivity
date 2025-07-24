
% create leadfields for spinal cord

%% Load geometries

addpath('D:\torso_tools')

sub='OP00219';
analysis='merged';
save_dir = fullfile('D:\MSST001', [sub '_' analysis]);

geoms_path = fullfile('D:\MSST001', [sub '_' analysis]);
geoms = load(fullfile(geoms_path, 'geoms.mat'));

filename = fullfile('D:\MSST001', ...
    ['sub-' sub], ...
    'ses-001', ...
    'meg', ...
    ['pmergedoe1000mspddfflo45hi45hfcstatic_001_array1.mat']);
D=spm_eeg_load(filename);

ordering = {'wm','bone','heart','lungs','torso'};

tt_add_bem;
reduction_factor = 0.3;  

for ii = 1:numel(ordering)
    field = ['mesh_' ordering{ii}];
    mesh_tmp = geoms.(field);

    if ii <= 2
        patch_in.vertices = mesh_tmp.vertices;
        patch_in.faces = mesh_tmp.faces;

        patch_out = reducepatch(patch_in, reduction_factor);

        pos = patch_out.vertices;
        tri = patch_out.faces;
    else
        pos = mesh_tmp.vertices;
        tri = mesh_tmp.faces;
    end

    bnd(ii).pos = pos;
    bnd(ii).tri = tri;
    bnd(ii).unit = 'mm';

    orient = hbf_CheckTriangleOrientation(bnd(ii).pos, bnd(ii).tri);
    if orient == 2
        bnd(ii).tri = bnd(ii).tri(:, [1 3 2]);
    end

    bnd(ii) = ft_convert_units(bnd(ii), 'm');
end

%% Conductivities
cratio = 40;
ci = [0.33 0.33/cratio .62 .05 .23];
co = [0.33/cratio .23 .23 .23 0];

cfg = [];
cfg.method = 'bem_hbf';
cfg.conductivity = [ci;co];
vol = ft_prepare_headmodel(cfg, bnd);

%% Source and sensor setup
src = [];
src.pos = geoms.sources_center_line.pos;
src.inside = true(size(src.pos, 1), 1);
src.unit = 'mm';
src = ft_convert_units(src, 'm');

grad=D.sensors('MEG');
grad = ft_convert_units(grad, 'm');
%grad = ft_datatype_sens(grad);  

%% Compute leadfields (all sources and orientations in one call)
cfg = [];
cfg.grad = grad;
cfg.headmodel = vol;
cfg.sourcemodel = src;
cfg.reducerank = 'no';
cfg.channel = 'all';
cfg.normalize = 'no';
leadfield = ft_prepare_leadfield(cfg);

save(fullfile(save_dir,'leadfields'), "leadfield")

