%% Beamformer
%outputs time series for each of the 3 orientations and each sourcepoint in
%spinal cord and brain

%note that one subject havs filename that deviate from the template. 224
%use static_002

% subs = {'OP00212', 'OP00213', 'OP00214', 'OP00215', 'OP00219', ...
%     'OP00220', 'OP00221', 'OP00224', 'OP00225', 'OP00226'};


clear all
close all
clc

addpath('D:\brainspineconnectivity')
addpath('D:\spm')
spm('defaults','EEG')

subs={'OP00212'};
generic_dir = 'D:\new_leadfields_and_geom';

HFC=1;
which_ori='all';

for k=1:length(subs)

    sub=subs{k};

    save_dir = fullfile('D:\MSST001', [sub '_contrast']);

    if ~exist(save_dir, 'dir')
        mkdir(save_dir)
    end


    if HFC
        filename_merged = fullfile('D:\MSST001', ...
            ['sub-' sub], ...
            'ses-001', ...
            'meg', ...
            ['pmergedoe1000mspddfflo45hi45hfcstatic_001_array1.mat']);
    else
        filename_merged = fullfile('D:\MSST001', ...
            ['sub-' sub], ...
            'ses-001', ...
            'meg', ...
            ['pmergedoe1000mspddfflo45hi45static_001_array1.mat']);


    end


    % SETUP AND PARAMS
    lambda = 0.01;  % Regularization parameter
    covtime = [-Inf Inf];  % Covariance time window

    [Gx, Gy, Gz] = build_leadfield_matrices(fullfile(generic_dir,'cervical_complex_brainspine_sources'));
    %note that this outputs with chans in same order as grad struct!!
    %ie use grad labels

    n_sourcepoints=size(Gx,1);
    n_channels=size(Gx,2);

    geomfile = fullfile(generic_dir, 'geometries_realistic_bone_cervical.mat');
    castforward = load(geomfile);

    mtorso = castforward.mesh_torso;


    %%  Plot torso mesh

    figure('Color','w');
    ft_plot_mesh(mtorso, 'facecolor', [0.8 0.8 0.8], 'facealpha', 0.1, 'edgecolor', 'none');
    hold on
    camlight headlight; lighting gouraud; material dull

    % Plot brain sources
    plot3(castforward.sources_brain.pos(:,1), ...
        castforward.sources_brain.pos(:,2), ...
        castforward.sources_brain.pos(:,3), ...
        'o', 'MarkerSize', 2, 'MarkerFaceColor', [0 0.45 0.74], ...
        'MarkerEdgeColor', 'k', 'DisplayName','Brain sources');

    % Plot spinal sources
    plot3(castforward.sources_cent.pos(:,1), ...
        castforward.sources_cent.pos(:,2), ...
        castforward.sources_cent.pos(:,3), ...
        's', 'MarkerSize', 6, 'MarkerFaceColor', [0.85 0.33 0.1], ...
        'MarkerEdgeColor', 'k', 'DisplayName','Spinal sources');

    %% LOAD  DATA
    D = spm_eeg_load(filename_merged);

    % BEAMFORMER PROCESSING LOOP

    for nn = 1% (runs)
        Dbf = D;
        chanind = indchantype(Dbf, 'MEGMAG'); %all channels with pos (=grad)
        chanind = chanind(~ismember(chanind, Dbf.badchannels)); %exclude bad channels here
        channames = Dbf.chanlabels(chanind); %excl bad chans
        sensors = Dbf.sensors('MEG');

        chanpos_s = zeros(length(chanind), 3);

        for f = 1:length(chanind)
            idx = find(strcmp(sensors.label, channames{f}));
            chanpos_s(f,:) = sensors.coilpos(idx,:);
        end

        covsamples = find(Dbf.time > covtime(1) & Dbf.time < covtime(2));
        Ntrials = Dbf.ntrials;
        allDbf=Dbf(:,:,:);

        % Covariance matrix
        covdata = zeros(length(chanind));
        for f = 1:Ntrials
            data = squeeze(allDbf(chanind, covsamples, f));
            data = data - mean(data, 2);
            covdata = covdata + data * data';
        end
        covdata = covdata / Ntrials;
        invcovdata = pinv(covdata + lambda * trace(covdata) * eye(size(covdata)));

        % lead fields for brain and spinal cord
        lfind = find(ismember(sensors.label, Dbf.chanlabels(chanind)));
        weights = zeros(n_sourcepoints, length(chanind), 3);

        L=zeros(n_sourcepoints,n_channels,3);
        L(:,:,1)=Gx;
        L(:,:,2)=Gy;
        L(:,:,3)=Gz;

        for cpind = 1: n_sourcepoints
            H = squeeze(L(cpind, lfind,:));
            weights(cpind, :, :) = (pinv(H' * invcovdata * H) * H' * invcovdata)';
        end

        % Apply beamformer
        wdata = zeros( n_sourcepoints, Ntrials, length(covsamples), 3);
        tcord = zeros( n_sourcepoints, length(covsamples), 3);

        for cpind = 1:n_sourcepoints
            for f = 1:Ntrials
                data = squeeze(allDbf(chanind, covsamples, f));
                data = data - mean(data, 2);
                wdata(cpind, f, :, :) = data' * squeeze(weights(cpind, :, :));
            end

        end
        stidx=find(strcmp(D.conditions,'static'));
        reidx=find(strcmp(D.conditions,'rest'));

        wdatarest=wdata(:,reidx,:,:);
        wdatastatic=wdata(:,stidx,:,:);

        % Save results
        %         if HFC
        %             bffilename = fullfile(save_dir, ['bfdata_', sub]);
        %
        %         else
        %             bffilename = fullfile(save_dir, ['bfdata_nohfc_', sub]);
        %         end
        %         spmfilename = Dbf.fullfile;
        %         save(bffilename, 'wdata', 'weights', 'covsamples', 'lambda', ...
        %             'spmfilename',  'chanind', 'wdatarest', 'wdatastatic')%,'brainsrc','headmesh','source_brain');


    end

    dat=spm2fieldtrip(D);

    powstat=zeros(D.ntrials,n_sourcepoints); %save power for each source pt

    % choose orientation or sum over orientations
    if strcmp(which_ori, 'all')
        ori_list = 1:3;
    else
        ori_list = str2double(which_ori);
    end

    if strcmp(which_ori, 'all')
        powstat = 0;
        powrest = 0;
    end

    for o = ori_list


        data = [];
        data.fsample = dat.hdr.Fs;  % sampling frequency
        data.label   = arrayfun(@(x) sprintf('src%03d',x), 1:n_sourcepoints, 'uni',0);

        for j = 1:D.ntrials
            data.trial{j} = squeeze(wdata(:,j,:,o)); % [nSrc x nTime]
            data.time{j}  = (0:size(wdata,3)-1)./data.fsample;
        end


        % Spectral power analysis
        cfg = [];
        cfg.output     = 'pow';
        cfg.method     = 'mtmfft';
        cfg.foilim     = [15 30];
        cfg.tapsmofrq  = 2;
        cfg.keeptrials = 'yes';
        cfg.pad        = 'nextpow2';

        fdat = ft_freqanalysis(cfg, data);


        % Trim trials and compute mean power
        stattrials=find(strcmp(D.conditions,'static'));
        resttrials=find(strcmp(D.conditions,'rest'));

        ntrials = min([numel(stattrials) numel(resttrials)]);
        meanpowstat = mean(fdat.powspctrm(stattrials,:,:), 3);
        meanpowrest = mean(fdat.powspctrm(resttrials,:,:), 3);

        if strcmp(which_ori, 'all')
            powstat = powstat + meanpowstat;
            powrest = powrest + meanpowrest;
        else
            powstat = meanpowstat;
            powrest = meanpowrest;
        end
    end


    %% t tests and permutation testing

    n_permutations = 100;

    tvals = zeros(1, n_sourcepoints);
    pvals = zeros(1, n_sourcepoints);

    for k = 1:n_sourcepoints
        statdat = log(powstat(1:ntrials, k));
        restdat = log(powrest(1:ntrials, k));

        [~, p, ~, stats] = ttest2(statdat, restdat, 'Tail', 'left');
        warning('one-tailed dep test, testing for power decrease') %i.e. neg t vals
        tvals(k) = stats.tstat;
        pvals(k) = p;
    end

    null_tvals = zeros(n_sourcepoints, n_permutations);

    for i = 1:n_permutations
        for k = 1:n_sourcepoints
            combined = [log(powstat(1:ntrials, k)); log(powrest(1:ntrials, k))];
            shuffled = combined(randperm(length(combined)));

            % Split back into two groups
            group1 = shuffled(1:ntrials);
            group2 = shuffled(ntrials+1:end);

            % Perform t-test
            [~, ~, ~, stats] = ttest2(group1, group2, 'Tail', 'left');
            null_tvals(k, i) = stats.tstat;
        end
    end

    thresholds = prctile(null_tvals, 5, 2);  % lower 5% threshold (left tail)
    [pk, peakind] = min(tvals);

    %% visualise t-stats
    
    spinemesh=castforward.mesh_wm;
    src_spine=castforward.sources_cent;
    func_vals_spine=tvals(1:length(src_spine.pos));
    
    brainmesh=castforward.mesh_brain;
    src_brain=castforward.sources_brain;
    src_brain.unit='mm';
    func_vals_brain=tvals(length(src_spine.pos)+1:end);



source = [];
source.pos    = src_brain.pos;             % [nSources x 3] xyz positions
source.inside = true(length(src_brain.pos),1); 
source.pow    = func_vals_brain';              % functional values per source
source.thresh=thresholds(length(src_spine.pos)+1:end);

source_spine=[];
source_spine.pos=src_spine.pos;
source_spine.inside=true(length(src_spine.pos),1);
source_spine.pow=func_vals_spine';
source_spine.thresh=thresholds(1:length(src_spine.pos));

addpath('D:\brainspineconnectivity\plotting')

cfg = [];
cfg.parameter = {'pow', 'thresh'};
source_int = ft_sourceinterpolate(cfg, source, brainmesh); 
spine_int=ft_sourceinterpolate(cfg,source_spine, spinemesh);

 source_int.mask = (source_int.pow < 0) & (source_int.pow < source_int.thresh);
 spine_int.mask = (spine_int.pow < 0) & (spine_int.pow < spine_int.thresh);



%% source plots

minval = min([source_int.pow(source_int.mask));

figure
subplot(121)
cfg = [];
cfg.figure='gcf';
cfg.method = 'surface';
cfg.funparameter = 'pow';
% cfg.maskparameter = 'mask';
% cfg.maskstyle='opacity';
cfg.funcolormap = flipud(brewermap(512,'Blues')); 
cfg.funcolorlim = 'minzero';
cfg.projmethod = 'nearest';
cfg.surffile = brainmesh;   % your mesh struct with .pos and .tri
ft_sourceplot(cfg, source_int); 
view(180,-8)
camlight

subplot(122)
cfg2=cfg;
cfg.funcolorlim = 'minzero';
cfg.projmethod = 'nearest';
cfg2.surffile = spinemesh;   % your mesh struct with .pos and .tri
% cb2 = colorbar;
% cb2.Position = [0.85 0.15 0.02 0.3]; % manual placement
% ylabel(cb2,'t (spine)','FontSize',10)
ft_sourceplot(cfg2, spine_int);
view(80,25)
camlight


f=figure; plot(src_spine.pos(:,2),func_vals_spine, 'ko-');
hold on; plot(src_spine.pos(:,2), thresholds(1:length(src_spine.pos)),'k--','LineWidth',1)
box off
set(gca,'FontSize',14)
set(gca, 'YDir', 'reverse') %neg upwards
title(sprintf('subject %s',sub))



end

toc
%now we have time series for each source point (only good channels)

