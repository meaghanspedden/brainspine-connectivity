
function [D,fullfilename]=convert_opm2spm(datadir,posfile,sub,bids_session,task,MEGrun,bids_precision)


fprintf('\n Converting to SPM..\n for subject %s, session %s, task %s, run %s\n',sub,bids_session,task,MEGrun)

save_dir= [fullfile(datadir,'spm')];
mkdir(save_dir);
cd(save_dir);

%% Load the OPM data---------------------------------------------------

% load head data so we get separate grad structure
disp('Loading data...');
cfg             = [];
cfg.folder      = datadir;
cfg.precision   = bids_precision;
cfg.bids.task   = task;
cfg.bids.sub    = sub;
cfg.bids.ses    = bids_session;
cfg.bids.run    = MEGrun;
cfg.binpath     = datadir;
cfg.dsPrefix    = false;

S=[];
S.precision = bids_precision;


file_name_bids = ['sub-' cfg.bids.sub '_ses-' cfg.bids.ses ...
    '_task-' cfg.bids.task '_run-' cfg.bids.run '_meg.lvm'];
path_to_bin_file = fullfile(cfg.folder,['sub-' cfg.bids.sub],...
    ['ses-' cfg.bids.ses],'meg');
fullfilename=[path_to_bin_file filesep file_name_bids];




if isempty(posfile)
 
            warning('no positions file')
else
    S.positions=posfile;
end


S.data=[path_to_bin_file filesep file_name_bids];

if ~exist(S.data),
    error('Cannot find %s',fullfilename);
end;

D= spm_opm_create(S);



