
cfg                   = [];
cfg.covariance        = 'yes';
cfg.channel           = good_spchans;
tlock = ft_timelockanalysis(cfg, ftdat);

cfg                     = [];
cfg.method              = 'lcmv';
cfg.headmodel           = dummyvol;
cfg.sourcemodel.LF      = Lf.leadfield;
cfg.grid                = sources_cent;
cfg.unit                = Lf.unit;
cfg.lcmv.keepfilter     = 'yes';
cfg.channel             = Lf.label;
source_time = ft_sourceanalysis(cfg, tlock);

beamformer=source_time.avg.filter;

spine_source_data = [];
spine_source_data.label = arrayfun(@(x) sprintf('source%d', x), 1:nsourcepoints, 'UniformOutput', false);

spine_source_data.time = ftdat.time;

[is_in, idx_good] = ismember(good_spchans, ftdat.label);


for sp=1:nsourcepoints
    for i=1:length(ftdat.trial) %gives a time series for each ori
        spine_source_data.trial{i} = beamformer{sp} * ftdat.trial{i}(idx_good,:);
    end
end


timeseries = cat(2, spine_source_data.trial{:});
[u1, s1, v1] = svd(timeseries, 'econ');

virtualchanneldata = [];
virtualchanneldata.label = spine_source_data.label;
virtualchanneldata.time = ftdat.time;


    for k = 1:length(ftdat.trial)
        for sp=1:nsourcepoints
        virtualchanneldata.trial{k}(sp,:) = u1(:,1)' * beamformer{sp} * ftdat.trial{k}(idx_good,:);
        end
    end
savename=sprintf('VC_spine_%s',sub);
save(fullfile(save_dir,savename), 'virtualchanneldata')
