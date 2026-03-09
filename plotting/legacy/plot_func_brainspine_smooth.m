function plot_func_brainspine_smooth(mesh, sc_src, sc_func, br_src, br_func, method)
% Inputs:
% - mesh: surface mesh for anatomy
% - sc_src/br_src: source structs with .pos and .inside
% - sc_func/br_func: functional values for spinal cord / brain
% - method: 'smooth' or 'raw' (smooth for both brain and spinal cord if 'smooth')

% Colormap
colmaptouse = 'summer';%brewermap([], 'viridis');

% Get inside points
inside_sc = sc_src.pos;
% inside_br = br_src.pos(br_src.inside, :);

% Color limits
clim = [min(sc_func(:)) 0];%; br_func(:)]));

% Plot base mesh
figure;
ft_plot_mesh(mesh, 'FaceColor', [0.8 0.8 1], 'facealpha', 0.2, ...
    'EdgeColor', 'none', 'AmbientStrength', 0.15); hold on;

% Raw plots for both brain and spinal cord
ft_plot_cloud(inside_sc, sc_func, 'scalerad', 'no', ...
    'cloudtype', 'surf', 'radius', 6,  ...
    'colormap', colmaptouse); %'clim', clim

% ft_plot_cloud(br_src, br_func, 'scalerad', 'no', ...
%     'cloudtype', 'surf', 'radius', 4,...
% 'colormap', colmaptouse);

xlabel('X'); ylabel('Y'); zlabel('Z');
axis image;
lighting flat;           % Flat lighting for a matte look
material dull;           % Dull material reduces specular highlights
ylim([-119 213])
view(83,5)
camlight('headlight');   % Light follows the camera but less shiny than default camlight

colorbar
end
