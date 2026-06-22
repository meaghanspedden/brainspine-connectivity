%% sensitivity_sensor_shifts_RL.m
% Computes R² of leadfield correlation (RL dipole orientation) between
% original and each sensor-shifted geometry, across spinal cord positions.
% Visualises R² along the cord, grouped by shift magnitude band.

geoms_path = 'C:\Users\mspedden\Documents\for_meaghan\geom_sens_shift';
lf_path    = 'C:\Users\mspedden\Documents\for_meaghan\bem_fields';

% Bundle metadata
bundle_names = {'Small (~2 mm)', 'Medium (~5 mm)', 'Large (~10 mm)'};
bundle_cols  = [0.18 0.55 0.90;   % bundle 1 — blue
                0.98 0.60 0.20;   % bundle 2 — orange
                0.65 0.08 0.05];  % bundle 3 — red

%% Load original leadfield
orig_lf_dir  = fullfile(lf_path, 'geometries_sensor_original');
orig_lf_file = dir(fullfile(orig_lf_dir, 'leadfield_*.mat'));
assert(~isempty(orig_lf_file), 'Could not find original leadfield in %s', orig_lf_dir);

orig      = load(fullfile(orig_lf_dir, orig_lf_file(1).name));
LF_orig   = extractLF_RL(orig.leadfield_cord);
src_pos   = orig.leadfield_cord.pos(orig.leadfield_cord.inside, :);
cord_pos  = src_pos(:, 2) * 1000;   % inferior-superior axis, m -> mm
n_sources = size(LF_orig, 2);

%% Discover all sensor-shifted leadfield folders
lf_dirs = dir(fullfile(lf_path, 'geometries_sensor_bundle*'));
lf_dirs = lf_dirs([lf_dirs.isdir]);
fprintf('Found %d sensor-shifted conditions\n', numel(lf_dirs));

%% Loop over shifted conditions
results = struct();
row = 1;

for ci = 1:numel(lf_dirs)
    folder = lf_dirs(ci).name;
    label  = regexprep(folder, '^geometries_', '');

    lf_files = dir(fullfile(lf_path, folder, 'leadfield_*_bem_experimental.mat'));
    if isempty(lf_files)
        warning('No leadfield found in %s — skipping', folder); continue
    end

    geom_file = fullfile(geoms_path, [folder '.mat']);
    if ~exist(geom_file, 'file')
        warning('No geometry file found for %s — skipping', folder); continue
    end

    % Load shift metadata
    geom      = load(geom_file, 'shift_vec_mm');
    shift_vec = geom.shift_vec_mm;
    shift_mag = norm(shift_vec);

    % Parse bundle number from folder name
    tok = regexp(folder, 'bundle(\d+)', 'tokens', 'once');
    bundle_idx = str2double(tok{1});

    fprintf('Processing: %s  |shift|=%.2f mm (%d/%d)\n', ...
        label, shift_mag, ci, numel(lf_dirs));

    shifted  = load(fullfile(lf_path, folder, lf_files(1).name));
    LF_shift = extractLF_RL(shifted.leadfield_cord);

    % R² per source
    r2 = zeros(n_sources, 1);
    for k = 1:n_sources
        r = corrcoef(LF_orig(:,k), LF_shift(:,k));
        r2(k) = r(1,2)^2;
    end

    [min_r2, min_idx] = min(r2);
    below99           = find(r2 < 0.99, 1, 'first');
    below95           = find(r2 < 0.95, 1, 'first');

    results(row).label      = label;
    results(row).bundle_idx = bundle_idx;
    results(row).shift_vec  = shift_vec;
    results(row).shift_mag  = shift_mag;
    results(row).r2_vector  = r2;
    results(row).median_r2  = median(r2);
    results(row).min_r2     = min_r2;
    results(row).min_pos    = src_pos(min_idx, :);
    results(row).below99_idx = below99;
    results(row).below95_idx = below95;

    fprintf('  median R²=%.4f  min R²=%.4f\n', median(r2), min_r2);
    row = row + 1;
end

save(fullfile(lf_path, 'sensitivity_results_sensor_RL.mat'), ...
    'results', 'src_pos', 'cord_pos', 'bundle_names', 'bundle_cols');
fprintf('\nDone. %d conditions processed.\n', row-1);

%% ---- R² figure ----
figure('Color','w');
hold on;

for bi = 1:3
    bundle_res = results([results.bundle_idx] == bi);
    col        = bundle_cols(bi,:);
    for ki = 1:numel(bundle_res)
        if ki == 1
            plot(cord_pos, bundle_res(ki).r2_vector, ...
                'Color', [col 0.6], 'LineWidth', 1.5, ...
                'DisplayName', bundle_names{bi});
        else
            plot(cord_pos, bundle_res(ki).r2_vector, ...
                'Color', [col 0.6], 'LineWidth', 1.5, ...
                'HandleVisibility', 'off');
        end
    end
end

% Threshold lines
yline(0.99, 'k:', 'LineWidth', 1.0, 'HandleVisibility', 'off');
yline(0.95, 'k--', 'LineWidth', 1.0, 'HandleVisibility', 'off');
text(max(cord_pos), 0.990, ' R²=0.99', 'FontSize', 9, ...
    'Color', [0.3 0.3 0.3], 'VerticalAlignment', 'bottom');
text(max(cord_pos), 0.950, ' R²=0.95', 'FontSize', 9, ...
    'Color', [0.3 0.3 0.3], 'VerticalAlignment', 'bottom');

ylim([0.8 1.02]);
pad = 0.05 * range(cord_pos);
xlim([min(cord_pos)-pad  max(cord_pos)+pad]);
xlabel('Position along cord (mm, inferior → superior)');
ylabel('R²');
title('Sensor shift sensitivity — R/L dipole orientation', ...
    'FontWeight', 'normal', 'FontSize', 12);
legend('Location', 'southwest', 'FontSize', 10);
grid on; box on;

%% ---- Sensor shift visualisation ----
orig_geom = load(fullfile(geoms_path, 'geometries_sensor_original.mat'));
orig_sens  = orig_geom.experimental_sensors.chanpos;
torso      = orig_geom.mesh_torso;

geom_files = dir(fullfile(geoms_path, 'geometries_sensor_bundle*.mat'));

figure('Color','w','Position',[100 100 900 700]);
hold on; axis equal; grid on;

% Torso mesh
trisurf(torso.faces, ...
    torso.vertices(:,1), torso.vertices(:,2), torso.vertices(:,3), ...
    'FaceColor', [0.85 0.75 0.65], 'FaceAlpha', 0.10, ...
    'EdgeColor', [0.6 0.5 0.4],    'EdgeAlpha', 0.06, ...
    'HandleVisibility', 'off');

% Original sensors
scatter3(orig_sens(:,1), orig_sens(:,2), orig_sens(:,3), 35, ...
    'MarkerFaceColor', [0.2 0.2 0.2], 'MarkerEdgeColor', 'none', ...
    'DisplayName', 'Original');

for bi = 1:3
    col        = bundle_cols(bi,:);
    bfiles     = geom_files(contains({geom_files.name}, sprintf('bundle%d', bi)));
    bfiles     = sort_by_shift_idx(bfiles);
    n_shifts   = numel(bfiles);
    alphas     = linspace(0.3, 0.9, n_shifts);

    for si = 1:n_shifts
        g        = load(fullfile(geoms_path, bfiles(si).name));
        sens_pos = g.experimental_sensors.chanpos;

        scatter3(sens_pos(:,1), sens_pos(:,2), sens_pos(:,3), 15, ...
            'MarkerFaceColor', col, 'MarkerEdgeColor', 'none', ...
            'MarkerFaceAlpha', alphas(si), ...
            'HandleVisibility', 'off');

        % Arrows on subsample
%         K   = max(1, floor(size(orig_sens,1)/10));
%         idx = 1:K:size(orig_sens,1);
%         for k = idx
%             dp = sens_pos(k,:) - orig_sens(k,:);
%             if norm(dp) > 0.01
%                 quiver3(orig_sens(k,1), orig_sens(k,2), orig_sens(k,3), ...
%                     dp(1), dp(2), dp(3), 0, ...
%                     'Color', [col alphas(si)], 'LineWidth', 1.0, ...
%                     'HandleVisibility', 'off');
%             end
%         end
    end

    % One legend entry per bundle
    scatter3(nan, nan, nan, 30, ...
        'MarkerFaceColor', col, 'MarkerEdgeColor', 'none', ...
        'DisplayName', bundle_names{bi});
end

legend('Location', 'bestoutside', 'FontSize', 10);
xlabel('x (mm)'); ylabel('y (mm)'); zlabel('z (mm)');
title('Sensor array shifts by magnitude band', ...
    'FontWeight', 'normal', 'FontSize', 12);
view(3); lighting gouraud; camlight('headlight');
rotate3d on; hold off;

%% Helper — extract RL orientation (col 1) across all sensors
function LF = extractLF_RL(leadfield_cord)
    idx       = find(leadfield_cord.inside);
    n_sources = numel(idx);
    n_sensors = numel(leadfield_cord.label);
    LF        = zeros(n_sensors, n_sources);
    for k = 1:n_sources
        LF(:, k) = leadfield_cord.leadfield{idx(k)}(:, 1);
    end
end

%% Helper — sort struct array by shift index number in filename
function sorted = sort_by_shift_idx(file_struct)
    names = {file_struct.name};
    idxs  = zeros(1, numel(names));
    for k = 1:numel(names)
        tok      = regexp(names{k}, 'shift(\d+)', 'tokens', 'once');
        idxs(k)  = str2double(tok{1});
    end
    [~, order] = sort(idxs);
    sorted = file_struct(order);
end