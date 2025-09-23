function [Gx, Gy, Gz] = build_leadfield_matrices(folder_path, source_option)
    % BUILD_LEADFIELD_MATRICES formats leadfield matrices from saved .mat files
    %
    % Usage:
    %   [Gx, Gy, Gz] = build_leadfield_matrices(folder_path, source_option)
    %
    % Inputs:
    %   folder_path   - path containing spine_*_ori*.mat and brain_*_ori*.mat files
    %   source_option - 'spine' (only spinal sources) or 'both' (spine + brain)
    %
    % Outputs:
    %   Gx, Gy, Gz    - leadfield matrices for X, Y, Z orientations

    if nargin < 2
        source_option = 'both'; % default
    end

    % Collect files depending on option
    switch lower(source_option)
        case 'spine'
            files = dir(fullfile(folder_path, 'spine_*_ori*.mat'));
        case 'both'
            files = [dir(fullfile(folder_path, 'spine_*_ori*.mat')); ...
                     dir(fullfile(folder_path, 'brain_*_ori*.mat'))];
        otherwise
            error('source_option must be ''spine'' or ''both''.');
    end

    if isempty(files)
        error('No matching files found in the specified folder.');
    end

    % Parse filenames
    nF = numel(files);
    groups = cell(nF,1);   % 'spine' or 'brain'
    src_idx = zeros(nF,1);
    ori_idx = zeros(nF,1);

    for i = 1:nF
        tok = regexp(files(i).name, '^(spine|brain)_(\d+)_ori(\d+)\.mat$', 'tokens', 'once');
        if isempty(tok)
            error('Filename "%s" does not match expected format.', files(i).name);
        end
        groups{i} = tok{1};
        src_idx(i) = str2double(tok{2});  
        ori_idx(i) = str2double(tok{3});  % expected 1..3 for X/Y/Z
    end

    % Unique source lists per group
    spine_mask = strcmp(groups, 'spine');
    brain_mask = strcmp(groups, 'brain');
    spine_sources = unique(src_idx(spine_mask));
    brain_sources = unique(src_idx(brain_mask));

    % Sort them
    spine_sources = sort(spine_sources, 'ascend');
    brain_sources = sort(brain_sources, 'ascend');

    n_spine = numel(spine_sources);
    n_brain = numel(brain_sources);
    n_sourcepoints = n_spine + n_brain;

    % Determine channels from one file
    tmp = load(fullfile(folder_path, files(1).name), 'Ls_x');
    n_sensors = size(tmp.Ls_x, 1);
    n_channels = n_sensors * 3;

    % Initialize outputs
    Gx = zeros(n_sourcepoints, n_channels);
    Gy = zeros(n_sourcepoints, n_channels);
    Gz = zeros(n_sourcepoints, n_channels);

    % Build index maps
    for i = 1:nF
        fname = fullfile(folder_path, files(i).name);
        data = load(fname, 'Ls_x', 'Ls_y', 'Ls_z'); % channel orientations

        % Row = [all X, then all Y, then all Z]
        row_data = [data.Ls_x(:)', data.Ls_y(:)', data.Ls_z(:)'];

        % Map to combined index
        if spine_mask(i)
            row = find(spine_sources == src_idx(i), 1, 'first');
        else % brain
            row = n_spine + find(brain_sources == src_idx(i), 1, 'first');
        end

        % Place into the correct matrix by orientation
        switch ori_idx(i)
            case 1
                Gx(row, :) = row_data;
            case 2
                Gy(row, :) = row_data;
            case 3
                Gz(row, :) = row_data;
            otherwise
                error('Unexpected orientation index %d in file %s', ori_idx(i), files(i).name);
        end
    end
end