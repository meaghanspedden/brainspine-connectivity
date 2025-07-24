% Set the directory containing the .lvm files
folderPath = 'D:\MSST001\sub-OP00221\ses-001\meg';  % <-- Change this to your target folder

% Get list of .lvm files
fileList = dir(fullfile(folderPath, '*.lvm'));

% Loop through each file
for i = 1:length(fileList)
    oldName = fileList(i).name;

    % Match pattern like 'rest_run-004_array1.lvm' or 'static_run-005_array1.lvm'
    pattern = '^(rest|static)_run-(00\d+_array\d+)\.lvm$';
    tokens = regexp(oldName, pattern, 'tokens');

    if ~isempty(tokens)
        prefix = tokens{1}{1};      % 'rest' or 'static'
        suffix = tokens{1}{2};      % '004_array1' or similar

        newName = [prefix '_' suffix '.lvm'];

        % Rename file
        movefile(fullfile(folderPath, oldName), fullfile(folderPath, newName));
        fprintf('Renamed: %s -> %s\n', oldName, newName);
    end
end

