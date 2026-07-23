function paths = brainspine_config()
%% brainspine_config.m
% Machine-specific paths for external dependencies. Edit the values below
% to match your local installation before running any pipeline script.
% Repo-internal paths (source/, plotting/, preprocessing/, and their
% subfolders) are NOT set here — each script derives them at runtime from
% its own location, so they need no editing regardless of where the repo
% is cloned.

paths.fieldtrip_path  = 'C:\Users\mspedden\Documents\fieldtrip';
paths.spm_path        = 'C:\Users\mspedden\Documents\spm';
paths.neurospec_path  = 'C:\Users\mspedden\Documents\neurospec211NEW\neurospec211';

paths.data_root       = 'C:\spinecoh_data';
paths.leadfields_dir  = 'C:\Leadfields meshes';
paths.save_dir        = 'C:\Users\mspedden\Documents\brainspine_savetest';
paths.figures_out_dir = 'C:\forGeneratingFigures';   % GET_FIGURES.m output

paths.geomfile   = fullfile(paths.leadfields_dir, 'geometries_experimental_withbrain.mat');
paths.T_mat_path = fullfile(paths.leadfields_dir, 'T.mat');
paths.lf_path    = fullfile(paths.leadfields_dir, 'leadfield_experimental_bslaw_experimental.mat');

end
