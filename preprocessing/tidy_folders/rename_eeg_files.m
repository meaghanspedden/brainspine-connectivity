% Set the folder path
folderPath = 'D:\OP00224_experiment\EMG';

% Define extensions to process
extensions = {'.vhdr', '.vmrk', '.eeg'};

% Loop through each extension
for extIdx = 1:length(extensions)
    ext = extensions{extIdx};
    
    % Get list of files for this extension
    fileList = dir(fullfile(folderPath, ['*' ext]));

    for i = 1:length(fileList)
        oldName = fileList(i).name;

        % Match pattern: anything before "_", then 3 digits, then extension
        pattern = '^(.*_)(\d{3})(\.[^.]+)$';
        tokens = regexp(oldName, pattern, 'tokens');

        if ~isempty(tokens)
            % Extract filename parts
            prefix = tokens{1}{1};  % e.g., 'OP00224_'
            digits = tokens{1}{2};  % e.g., '012'
            suffix = tokens{1}{3};  % e.g., '.vhdr'

            % Construct new filename with '0' inserted
            newName = [prefix '0' digits suffix];

            % Rename the file
            movefile(fullfile(folderPath, oldName), fullfile(folderPath, newName));
            fprintf('Renamed: %s → %s\n', oldName, newName);
        else
            fprintf('Skipped (no match): %s\n', oldName);
        end
    end
end
