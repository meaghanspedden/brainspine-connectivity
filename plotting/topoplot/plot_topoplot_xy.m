function plot_topoplot_xy(ax, positions, leadfield_data, cmap, min_val, max_val)
    % Interpolate sensor data onto a grid for smoother visualization (x-y plane)
    grid_res = 100;  % Resolution of the grid
    [grid_x, grid_y] = meshgrid( ...
        linspace(min(positions(:, 1)), max(positions(:, 1)), grid_res), ...
        linspace(min(positions(:, 2)), max(positions(:, 2)), grid_res) ...
    );

    interpolated_values = griddata( ...
        positions(:, 1), positions(:, 2), leadfield_data, ...
        grid_x, grid_y, 'natural' ...
    );

    % Plot topoplot using contourf
    contourf(ax, grid_x, grid_y, interpolated_values, 30, 'LineStyle', 'none');
    colormap(ax, cmap);

    % Use consistent limits for all tiles
    caxis(ax, [min_val, max_val]);

    % Add colorbar with unit label
    cb = colorbar(ax);
    cb.Label.String = 'fT';  % Label the colorbar with "fT"

    axis(ax, 'equal');
    axis(ax, 'off');
end
