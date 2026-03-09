function plot_func_brainspine_interp(mesh, sc_src, sc_func, br_src, br_func, method)
    % Inputs:
    % - mesh: surface mesh for anatomy
    % - sc_src/br_src: source structures with .pos and .inside
    % - sc_func/br_func: functional values for spinal cord / brain
    % - method: 'interp' or 'raw' (only brain uses interpolation here)

    % Colormap
    colmaptouse ='parula';% brewermap([], 'PrGn');
    
    % Get inside points
    inside_sc = sc_src.pos(sc_src.inside, :);
    inside_br = br_src.pos(br_src.inside, :);

    % Color limits (symmetric)
    clim = max(abs([sc_func(:); br_func(:)]));

    % Plot base mesh
    figure;
    ft_plot_mesh(mesh, 'FaceColor', [0.8 0.8 1], 'facealpha', 0.2, ...
        'EdgeColor', 'none', 'AmbientStrength', 0.15); hold on;

    % Always plot spinal cord raw point cloud
    ft_plot_cloud(inside_sc, sc_func, 'scalerad', 'no', ...
        'cloudtype', 'surf', 'radius', 5, 'clim', [-clim clim], ...
        'colormap', flipud(colormap(colmaptouse)));

    if strcmp(method, 'interp')
        % Interpolate brain only

        % Create interpolation grid based on brain points
        [Xq, Yq, Zq] = ndgrid( ...
            linspace(min(inside_br(:,1)), max(inside_br(:,1)), 50), ...
            linspace(min(inside_br(:,2)), max(inside_br(:,2)), 50), ...
            linspace(min(inside_br(:,3)), max(inside_br(:,3)), 50));

        % Interpolate brain functional data
        try
            Fbr = scatteredInterpolant(inside_br(:,1), inside_br(:,2), inside_br(:,3), br_func', 'natural', 'none');
        catch
            warning('Brain interpolation failed; falling back to linear method.');
            Fbr = scatteredInterpolant(inside_br(:,1), inside_br(:,2), inside_br(:,3), br_func', 'linear', 'none');
        end

        Vq_br = Fbr(Xq, Yq, Zq);

        % Plot brain isosurface
        p2 = patch(isosurface(Xq, Yq, Zq, Vq_br, 0.1 * clim));
        set(p2, 'FaceColor', [0.8 0.2 0.2], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
    else
        % Plot brain raw points if not interpolating
        ft_plot_cloud(inside_br, br_func, 'scalerad', 'no', ...
            'cloudtype', 'surf', 'radius', 4, 'clim', [-clim clim], ...
            'colormap', flipud(colormap(colmaptouse)));
    end

    xlabel('X'); ylabel('Y'); zlabel('Z');
    axis image; colorbar;
    camlight; lighting gouraud;
end
