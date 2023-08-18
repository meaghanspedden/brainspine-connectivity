function dataOut = getMetaData (subID, mydir)

if strcmp(subID,'116')
    metadat.filenames=strvcat([mydir '\sub-OP00116\ses-001\meg\pprhandoe1000msfdfflo45hi5sub-OP00116_ses-001_task-staticright_run-001_meg.mat'],...
    [mydir '\sub-OP00116\ses-001\meg\pplhandoe1000msfdfflo45hi5sub-OP00116_ses-001_task-staticleft_run-001_meg.mat']);

    metadat.posfile=[mydir '\sub-OP00116\ses-001\meg\sub-OP00116_ses-001_positions_new.tsv'];
    metadat.backshape=['D:\OP00116_experiment\cast space\OP0016_seated.stl'];
    metadat.cyl=['D:\OP00116_experiment\cast space\cylinder_good_space.stl'];
    metadat.brainchan_labels={'19','1B', 'MW', 'MY','OK','MZ','35','DI','MU','DQ','OH','DG','A1'};
    metadat.badchans={};
    metadat.dpath=[mydir '\sub-OP00116\ses-001\meg'];
    metadat.savepath=[mydir '\Coh_results00116'];
    metadat.layoutfile='D:\OP00116_experiment\lay_head_116';


elseif strcmp(subID,'122')
    metadat.filenames=strvcat([mydir '\sub-OP00122\ses-001\meg\pprhandoe1000msfdfflo45hi5ds_sub-OP00122_ses-001_task-static_right_run-002_meg.mat'],...
        [mydir '\sub-OP00122\ses-001\meg\pplhandoe1000msfdfflo45hi5ds_sub-OP00122_ses-001_task-static_left_run-001_meg.mat']);
    metadat.posfile=[mydir '\sub-OP00122\ses-001\meg\ds_sub-OP00122_ses-001_positions.tsv'];
    metadat.backshape=['D:\OP00122_experiment\cast space\OP0015_seated.stl'];
    metadat.cyl=['D:\OP00122_experiment\cast space\cylinder_good_space.stl'];
    metadat.brainchan_labels={'19','DG', 'OH', 'A1','1B', 'A9','JS','A6','DJ','MY','DS','OK','MI','17'};
    metadat.badchans={'ML-X','ML-Y','ML-Z','K4-Z','K4-X'};
    metadat.dpath=[mydir '\sub-OP00122\ses-001\meg'];
    metadat.savepath=[mydir '\Coh_results00122'];
    metadat.layoutfile='D:\OP00122_experiment\lay_head_122';



elseif strcmp(subID,'123')
    metadat.filenames=strvcat([mydir 'sub-OP00124\ses-001\meg\pprhandoe1000msfdfflo45hi5ds_sub-OP00124_ses-001_task-staticright_run-001_meg'],...
        [mydir '\sub-OP00124\ses-001\meg\pplhandoe1000msfdfflo45hi5ds_sub-OP00124_ses-001_task-staticleft_run-001_meg']);
    metadat.posfile=[mydir 'sub-OP00124\ses-001\meg\ds_sub-OP00124_ses-001_positions.tsv'];
    metadat.backshape=['D:\OP00123_experiment\cast space\OP00123.stl'];
    metadat.cyl=['D:\OP00123_experiment\cast space\cylinder_good_space.stl'];
    metadat.brainchan_labels={'19','DG', 'OH', 'A1','1B', 'A9','JS','A6','DJ','MY','DS','OK','MI','17'};
    metadat.badchans={'ML-X','ML-Y','ML-Z','K4-Z','K4-X'};
    metadat.dpath=[mydir '\sub-OP00124\ses-001\meg'];
    metadat.savepath=[mydir '\Coh_results00123'];
    metadat.layoutfile='D:\OP00123_experiment\lay_head_123';


else
    error('no subject with this ID')
end

dataOut=metadat;

end %function