%% format lead fields

function [Gx, Gy, Gz] = build_leadfield_matrices(folder_path)
    % Get list of all spine_*_ori*.mat files in the folder
    files = dir(fullfile(folder_path, 'spine_*_ori*.mat'));
    
    if isempty(files)
        error('No matching files found in the specified folder.');
    end
    
    % Parse source indices and orientations from filenames
    source_indices = [];
    ori_indices = [];
    for i = 1:length(files)
        name_parts = regexp(files(i).name, 'spine_(\d+)_ori(\d+)\.mat', 'tokens');
        if isempty(name_parts)
            error('Filename "%s" does not match expected format.', files(i).name);
        end
        source_indices(end+1) = str2double(name_parts{1}{1}); %#ok<AGROW>
        ori_indices(end+1)    = str2double(name_parts{1}{2}); %#ok<AGROW>
    end
    
    % Determine number of sources and sensors
    n_sourcepoints = max(source_indices);
    
    % Load one file to determine number of sensors
    tmp = load(fullfile(folder_path, files(1).name), 'Ls_x');
    n_sensors = size(tmp.Ls_x, 1);
    n_channels = n_sensors * 3;
    
    % Initialise output matrices
    Gx = zeros(n_sourcepoints, n_channels);
    Gy = zeros(n_sourcepoints, n_channels);
    Gz = zeros(n_sourcepoints, n_channels);
    
    % Loop through files and fill matrices
    for i = 1:length(files)
        fname = fullfile(folder_path, files(i).name);
        data = load(fname, 'Ls_x', 'Ls_y', 'Ls_z');
        
        % Build ordered row: all X, then all Y, then all Z
        row_data = [data.Ls_x(:)', data.Ls_y(:)', data.Ls_z(:)'];
        
        % Identify source point and orientation
        name_parts = regexp(files(i).name, 'spine_(\d+)_ori(\d+)\.mat', 'tokens');
        src_idx = str2double(name_parts{1}{1});
        ori_idx = str2double(name_parts{1}{2});
        
        % Assign to correct matrix
        switch ori_idx
            case 1
                Gx(src_idx, :) = row_data;
            case 2
                Gy(src_idx, :) = row_data;
            case 3
                Gz(src_idx, :) = row_data;
            otherwise
                error('Unexpected orientation index %d in file %s', ori_idx, files(i).name);
        end
    end
end