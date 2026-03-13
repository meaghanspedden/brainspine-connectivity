folderPath = 'D:\VCG\sub-OP00202\ses-001\meg'; % Change this to your directory
cd(folderPath)
files = dir(fullfile(folderPath, 'MNS-run-*_array*.*'));
subID = 'OP00202'; % Change to your subject ID
error('fix channels tsv?')
for i = 1:length(files)
    oldName = files(i).name;
    oldPath = fullfile(folderPath, oldName);
    
    % Extract run number
    runMatch = regexp(oldName, 'MNS-run-(\d+)_array', 'tokens');
    if isempty(runMatch)
        continue;
    end
    runNum = str2double(runMatch{1}{1});
    newRunStr = sprintf('%03d', runNum); % Ensure 3-digit format
    
    % Extract file extension
    [~, ~, ext] = fileparts(oldName);
    
    % Construct new filename
    newName = sprintf('sub-%s_ses-001_task-mns_run-%s_meg%s', subID, newRunStr, ext);
    newPath = fullfile(folderPath, newName);
    
    % Rename file
    movefile(oldPath, newPath);
    fprintf('Renamed: %s -> %s\n', oldName, newName);
end
