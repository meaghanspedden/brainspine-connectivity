


folderPath='D:\brain_spinalcord_opm_data\sub-OP00226\ses-001\meg';
% Regex pattern for allowed filenames (no extension in this match)
keepPattern = '^(static|rest)_00[1-9]_array1$';

% Get all items in the folder
items = dir(folderPath);

for i = 1:numel(items)
    if items(i).isdir
        % Skip '.' and '..', otherwise delete folder
        if ~ismember(items(i).name, {'.','..'})
            folderToDelete = fullfile(folderPath, items(i).name);
            rmdir(folderToDelete, 's'); % 's' removes contents recursively
            fprintf('Deleted folder: %s\n', folderToDelete);
        end
    else
        % Process files
        [~, name, ext] = fileparts(items(i).name);

        % Check keep conditions
        isKeepName = ~isempty(regexp(name, keepPattern, 'once'));
        isKeepTSV  = strcmpi(ext, '.tsv');

        % Delete if not matching keep criteria
        if ~(isKeepName || isKeepTSV)
            delete(fullfile(folderPath, items(i).name));
            fprintf('Deleted file: %s\n', items(i).name);
        end
    end
end

disp('Folder cleanup complete.');