% calc RMS error for each subject

sub123_targets=[63 62];
sub122_targets=[58 62];
sub116_targets=[67 52];


rightfiles={'07', '09','11','13'};
leftfiles={'23', '25', '27'};

EMGpath='D:\OP00116_experiment\EMGfiles';
EMGfiletemplate='000116_static_';
EMGfilename=[EMGfiletemplate,cell2mat(rightfiles(1)),'.smrx'];
EMGchanname='EMG 1'

addpath(genpath('D:\brainspineconnectivity'))


[~, EMGstruct]= readSpikeData(fullfile(EMGpath, EMGfilename),EMGchanname);
