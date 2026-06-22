%% CREATE SHIFTED SOURCE MODELS FOR SENSITIVITY ANALYSIS
% Generates shifted versions of the spinal cord source model along three
% anatomically meaningful axes:
%   - Left-Right:      global X axis
%   - Ventral-Dorsal:  global Y axis
%   - Rostral-Caudal:  local tangent along cord
%
% Shift magnitudes: +/-2, +/-4, +/-6 mm
% Total: 3 directions x 3 magnitudes x 2 signs = 18 shifted geometries

%% Compute local tangent along cord
pos = spine_sources.pos;   % Nx3, assumed ordered along cord
N   = size(pos, 1);

tangent = zeros(N, 3);
tangent(2:end-1,:) = pos(3:end,:) - pos(1:end-2,:);
tangent(1,:)       = pos(2,:)     - pos(1,:);
tangent(end,:)     = pos(end,:)   - pos(end-1,:);
tangent            = tangent ./ vecnorm(tangent, 2, 2);

%% Define shift conditions
mags = [2, 4, 6];
signs = [1, -1];
sign_labels = {'pos', 'neg'};

directions = struct();
directions(1).name = 'LR';
directions(1).type = 'cardinal';
directions(1).vec  = repmat([1 0 0], N, 1);

directions(2).name = 'VD';
directions(2).type = 'cardinal';
directions(2).vec  = repmat([0 1 0], N, 1);

directions(3).name = 'RC';
directions(3).type = 'local';
directions(3).vec  = tangent;   % Nx3 local tangent at each source point

%% Save original geometry
geom_original = struct();
geom_original.mesh_wm             = mesh_wm;
geom_original.mesh_bone           = mesh_bone;
geom_original.mesh_heart          = mesh_heart;
geom_original.mesh_lungs          = mesh_lungs;
geom_original.mesh_torso          = mesh_torso;
geom_original.sources_cent        = spine_sources;
geom_original.experimental_sensors = exp_sensors;

outfile_original = fullfile(savepath, 'geometries_original.mat');
save(outfile_original, '-struct', 'geom_original', '-v7.3');
fprintf('Saved: geometries_original.mat\n');

%% Generate and save each shifted geometry
for di = 1:numel(directions)
    vec  = directions(di).vec;    % Nx3
    name = directions(di).name;

    for mi = 1:numel(mags)
        for si = 1:numel(signs)

            shift_mm    = mags(mi) * signs(si);
            sign_label  = sign_labels{si};
            label       = sprintf('%s_%s%dmm', name, sign_label, mags(mi));

            fprintf('Creating shifted geometry: %s...\n', label);

            % Apply shift — for cardinal directions this is a uniform
            % translation; for RC it follows the local cord tangent
            shifted_sources     = spine_sources;
            shifted_sources.pos = pos + shift_mm * vec;

            % Package geometry — meshes and sensors unchanged
            geom_shifted = struct();
            geom_shifted.mesh_wm              = mesh_wm;
            geom_shifted.mesh_bone            = mesh_bone;
            geom_shifted.mesh_heart           = mesh_heart;
            geom_shifted.mesh_lungs           = mesh_lungs;
            geom_shifted.mesh_torso           = mesh_torso;
            geom_shifted.sources_cent         = shifted_sources;
            geom_shifted.experimental_sensors = exp_sensors;

            outfile = fullfile(savepath, ['geometries_' label '.mat']);
            save(outfile, '-struct', 'geom_shifted', '-v7.3');
            fprintf('  Saved: geometries_%s.mat\n', label);

        end
    end
end

fprintf('\nDone. Saved %d shifted geometries + 1 original.\n', ...
    numel(directions) * numel(mags) * numel(signs));