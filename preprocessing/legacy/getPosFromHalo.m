%% convert results from halo to positions tsv for spm opm create


folder = 'D:\VCG\sub-OP00203\ses-001\meg\Halo-run-001_19-02-2025_11-39-35\RESULTS FOR Halo-run-001_array1\Results 2025-02-20 11-06';
file = 'RESULTS_Localization Solver-2025-02-20 11-06.csv';
fullFilePath = fullfile(folder, file);
opts = detectImportOptions(fullFilePath, 'FileType', 'text');
haloPos = readtable(fullFilePath, opts);
 
cfg = struct;
cfg.haloResults = fullFilePath;
cfg.output = 'D:\VCG\sub-OP00203\positions';                      
cfg.align = true;                                                
%cfg.lvm= 'D:\VCG\sub-OP00202\ses-001\meg\sub-OP00202_ses-001_task-mns_run-001_meg.lvm';
cfg.plot = true;                                                  
%cfg.threshold = 20;                                               
 
[positions, goodPositions] = positionsFromHalo(cfg);