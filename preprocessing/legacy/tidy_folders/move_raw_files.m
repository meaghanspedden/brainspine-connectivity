% move_lvm_files.m
% =========================================================================
% Recursively move all .lvm files from any sub-folder containing
% "rest" or "static" into the target D:\MSST001\sub-OP00221\ses-001\meg
% =========================================================================

% === User‐configurable paths ===
baseDir   = 'D:\MSST001\sub-OP00226\ses-001\meg';      % top of subject folder
targetDir = baseDir;%fullfile(baseDir, 'ses-001', 'meg');  % existing destination

% check destination exists
if ~exist(targetDir,'dir')
    error('Destination folder does not exist:\n%s', targetDir);
end

% === find all .lvm files under baseDir recursively ===
% Requires R2016b or later for the "**" wildcard
allLVM = dir(fullfile(baseDir, '**', '*.lvm'));

% === loop through and move matching files ===
moved = 0;
for ii = 1:numel(allLVM)
    fld = allLVM(ii).folder;          % full path to the file's folder
    [~, subname] = fileparts(fld);    % last folder name

    % check for "rest" or "static" in the folder name
    if contains(subname, 'rest',   'IgnoreCase', true) || ...
       contains(subname, 'static', 'IgnoreCase', true)

        src  = fullfile(fld, allLVM(ii).name);
        dst  = fullfile(targetDir, allLVM(ii).name);

        % move (will overwrite if same-named file exists)
        movefile(src, dst);
        fprintf('Moved: %s -> %s\n', src, dst);
        moved = moved + 1;
    end
end

fprintf('\nDone. Total files moved: %d\n', moved);
