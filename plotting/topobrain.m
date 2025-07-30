function topobrain(brainlay, sensor_values)


% Get X and Z positions
chanXZ = brainlay.pos;           % X-Z positions
zLabels = brainlay.label;                   % Labels

% Create interpolation grid
x_range = linspace(min(chanXZ(:,1)), max(chanXZ(:,1)), 200);
z_range = linspace(min(chanXZ(:,2)), max(chanXZ(:,2)), 200);
[Xgrid, Zgrid] = meshgrid(x_range, z_range);

% Interpolate sensor values
F = scatteredInterpolant(chanXZ(:,1), chanXZ(:,2), sensor_values, 'natural', 'none');
Vgrid = F(Xgrid, Zgrid);

% Plot interpolated image
%figure;
imagesc(x_range, z_range, Vgrid);
set(gca, 'YDir', 'normal');
colormap jet;
colorbar;
caxis([0 (max(max(Vgrid)))])
hold on;

% Overlay Z-channel positions
scatter(chanXZ(:,1), chanXZ(:,2), 40, 'k', 'filled');
text(chanXZ(:,1), chanXZ(:,2), zLabels, ...
    'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', 'FontSize', 7);

head_outline=brainlay.mask{1};
smoothed_outline = smooth_head_outline_spline(head_outline, 50, .5);
% Optional: plot head outline if provided
%if nargin > 2 && ~isempty(head_outline)
    plot(smoothed_outline(:,1), smoothed_outline(:,2), 'k-', 'LineWidth', 2);
%end

xlabel('X (mm)');
ylabel('Z (mm)');
axis image;


