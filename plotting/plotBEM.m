%% PLOT_BEM_GEOMETRY.m
% Methods figure showing torso surface, vertebrae, spinal cord (mesh_wm),
% and OPM sensor positions. Heart, lungs, brain and outer skull removed.

clear all
close all
clc

%% =========================================================================
%  USER CONFIG
%% =========================================================================

cfg.fieldtrip_path = 'C:\Users\mspedden\Documents\fieldtrip';
cfg.geomfile       = 'C:\Leadfields meshes\geometries_cervical_realistic.mat';
cfg.fig_dir        = 'C:\Users\mspedden\Documents\brainspine_save_bemv2\figures';
cfg.saveFig        = 1;

%% =========================================================================
%  SETUP
%% =========================================================================

addpath(cfg.fieldtrip_path)
ft_defaults

%% =========================================================================
%  LOAD GEOMETRY
%% =========================================================================

fprintf('Loading geometry file...\n')
geom = load(cfg.geomfile);

geom_vars = fieldnames(geom);
fprintf('Variables in geometry file:\n');
for v = 1:length(geom_vars)
    vn  = geom_vars{v};
    val = geom.(vn);
    if isstruct(val) && (isfield(val,'pos') || isfield(val,'vertices'))
        nv = 0; nf = 0;
        if isfield(val,'vertices'), nv = size(val.vertices,1); elseif isfield(val,'pos'), nv = size(val.pos,1); end
        if isfield(val,'faces'),    nf = size(val.faces,1);    elseif isfield(val,'tri'), nf = size(val.tri,1); end
        fprintf('  %s: mesh (%d vertices, %d faces)\n', vn, nv, nf);
    else
        fprintf('  %s: %s\n', vn, class(val));
    end
end

mesh_torso = convert_mesh(geom.mesh_torso);
mesh_bone  = convert_mesh(geom.mesh_bone);
mesh_cord  = convert_mesh(geom.mesh_wm);

grad_mm    = geom.coils_3axis;

brain_idx  = grad_mm.chanpos(:,2) > 200;
spine_idx  = ~brain_idx;
brain_pos  = grad_mm.chanpos(brain_idx,:);
spine_pos  = grad_mm.chanpos(spine_idx,:);
fprintf('  Brain sensors: %d\n', sum(brain_idx))
fprintf('  Spine sensors: %d\n', sum(spine_idx))

%% =========================================================================
%  APPEARANCE
%% =========================================================================

% Torso: very transparent so cord and vertebrae are visible through it
% Bone:  semi-transparent so cord shows through
% Cord:  opaque so it reads clearly as the structure of interest
alpha_torso = 0.12;
alpha_bone  = 0.50;
alpha_cord  = 0.90;

col_torso  = [0.65 0.50 0.75];   % light purple
col_bone   = [0.85 0.80 0.95];   % pale lavender
col_cord   = [0.75 0.20 0.65];   % purple-pink — stands out against bone/torso

col_sensor = [0 0 0];
sensor_sz  = 25;

%% =========================================================================
%  FIGURE 1: Lateral view
%% =========================================================================

hfig_lat = figure('Color','w','Position',[100 100 600 700]);
hold on; axis equal; axis off;
view(90, 0);

ft_plot_mesh(mesh_torso, 'facecolor', col_torso, 'facealpha', alpha_torso, 'edgecolor', 'none');
ft_plot_mesh(mesh_bone,  'facecolor', col_bone,  'facealpha', alpha_bone,  'edgecolor', 'none');
ft_plot_mesh(mesh_cord,  'facecolor', col_cord,  'facealpha', alpha_cord,  'edgecolor', 'none');

scatter3(spine_pos(:,1), spine_pos(:,2), spine_pos(:,3), sensor_sz, col_sensor, 'filled');
scatter3(brain_pos(:,1), brain_pos(:,2), brain_pos(:,3), sensor_sz, col_sensor, 'filled');

camlight('headlight'); lighting gouraud;
title('Lateral view', 'FontSize', 14, 'FontWeight', 'normal');

h_torso = patch('Faces',[1 2 3],'Vertices',[0 0 0;1 0 0;0 1 0], 'FaceColor', col_torso, 'FaceAlpha', alpha_torso, 'EdgeColor', 'none');
h_bone  = patch('Faces',[1 2 3],'Vertices',[0 0 0;1 0 0;0 1 0], 'FaceColor', col_bone,  'FaceAlpha', alpha_bone,  'EdgeColor', 'none');
h_cord  = patch('Faces',[1 2 3],'Vertices',[0 0 0;1 0 0;0 1 0], 'FaceColor', col_cord,  'FaceAlpha', alpha_cord,  'EdgeColor', 'none');
h_sens  = scatter3(nan, nan, nan, sensor_sz, col_sensor, 'filled');
legend([h_torso h_bone h_cord h_sens], ...
    {'Torso', 'Vertebrae', 'Spinal cord', 'OPM sensors'}, ...
    'Location', 'southoutside', 'Orientation', 'horizontal', 'FontSize', 11, 'Box', 'off');

%% =========================================================================
%  FIGURE 2: Superior view
%% =========================================================================

hfig_sup = figure('Color','w','Position',[100 100 600 600]);
hold on; axis equal; axis off;
view(180, 90);

% Keep torso semi-transparent in superior view too so cord is visible
ft_plot_mesh(mesh_torso, 'facecolor', col_torso, 'facealpha', 0.30,      'edgecolor', 'none');
ft_plot_mesh(mesh_bone,  'facecolor', col_bone,  'facealpha', alpha_bone, 'edgecolor', 'none');
ft_plot_mesh(mesh_cord,  'facecolor', col_cord,  'facealpha', alpha_cord, 'edgecolor', 'none');

scatter3(spine_pos(:,1), spine_pos(:,2), spine_pos(:,3), sensor_sz, col_sensor, 'filled');
scatter3(brain_pos(:,1), brain_pos(:,2), brain_pos(:,3), sensor_sz, col_sensor, 'filled');

camlight('headlight'); lighting gouraud;
title('Superior view (sensor coverage)', 'FontSize', 14, 'FontWeight', 'normal');

h_torso2 = patch('Faces',[1 2 3],'Vertices',[0 0 0;1 0 0;0 1 0], 'FaceColor', col_torso, 'FaceAlpha', 0.30,      'EdgeColor', 'none');
h_bone2  = patch('Faces',[1 2 3],'Vertices',[0 0 0;1 0 0;0 1 0], 'FaceColor', col_bone,  'FaceAlpha', alpha_bone, 'EdgeColor', 'none');
h_cord2  = patch('Faces',[1 2 3],'Vertices',[0 0 0;1 0 0;0 1 0], 'FaceColor', col_cord,  'FaceAlpha', alpha_cord, 'EdgeColor', 'none');
h_sens2  = scatter3(nan, nan, nan, sensor_sz, col_sensor, 'filled');
legend([h_torso2 h_bone2 h_cord2 h_sens2], ...
    {'Torso', 'Vertebrae', 'Spinal cord', 'OPM sensors'}, ...
    'Location', 'southoutside', 'Orientation', 'horizontal', 'FontSize', 11, 'Box', 'off');

%% =========================================================================
%  SAVE
%% =========================================================================

if cfg.saveFig
    if ~exist(cfg.fig_dir,'dir'), mkdir(cfg.fig_dir); end
    savefig(hfig_lat, fullfile(cfg.fig_dir, 'BEM_geometry_lateral.fig'));
    print(hfig_lat,   fullfile(cfg.fig_dir, 'BEM_geometry_lateral.png'),  '-dpng', '-r300');
    savefig(hfig_sup, fullfile(cfg.fig_dir, 'BEM_geometry_superior.fig'));
    print(hfig_sup,   fullfile(cfg.fig_dir, 'BEM_geometry_superior.png'), '-dpng', '-r300');
    fprintf('Figures saved to %s\n', cfg.fig_dir);
end

fprintf('\nDone.\n')

%% =========================================================================
%  HELPER
%% =========================================================================

function mesh_out = convert_mesh(mesh_in)
    mesh_out = mesh_in;
    if isfield(mesh_in,'vertices') && ~isfield(mesh_in,'pos')
        mesh_out.pos = mesh_in.vertices;
    end
    if isfield(mesh_in,'faces') && ~isfield(mesh_in,'tri')
        mesh_out.tri = mesh_in.faces;
    end
end