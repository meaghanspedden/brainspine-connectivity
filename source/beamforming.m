%% Beamformer
%outputs time series for each of the 3 orientations and each sourcepoint in
%spinal cord

%note that two subjects have filenames that deviate from the template. 224
%use static_002 and 214 use 'spd' instead of 'spdd' (no ref) - so can't
%loop through all.

% subs = {'OP00212', 'OP00213', 'OP00214', 'OP00215', 'OP00219', ...
%     'OP00220', 'OP00221', 'OP00224', 'OP00225', 'OP00226'};


clear all
close all
clc
subs={'OP00214'};

generic_dir = 'D:\MSST001\generic_merged';

HFC=0;

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
    hstep = 10;  % Grid spacing for mesh sampling brain

    rootlf   = fullfile(generic_dir, 'leadfields.mat');
    geomfile = fullfile(generic_dir, 'geoms.mat');
    castforward = load(geomfile);


    % PLOT MESH & SOURCE POSITIONS
    mtorso = castforward.mesh_torso;
    spos = castforward.sources_center_line.pos;
    Ncordpoints = length(spos);


    % LOAD  DATA
    D = spm_eeg_load(filename_merged);
    ftdat=spm2fieldtrip(D);
    grad=D.sensors('MEG');
    labs=D.chanlabels;


    %% spinal cord---------------------------------------------------------------------

    load(rootlf) %leadfield
    num_sources=numel(leadfield.leadfield);
    num_channels=numel(leadfield.label);

    Gx = zeros(num_sources, num_channels);
    Gy = zeros(num_sources, num_channels);
    Gz = zeros(num_sources, num_channels);


    % Loop through each source
    for i = 1:num_sources
        L = leadfield.leadfield{i};  % n channels x 3 matrix
        Gx(i, :) = L(:, 1)';  % x-orientation
        Gy(i, :) = L(:, 2)';  % y-orientation
        Gz(i, :) = L(:, 3)';  % z-orientation
    end

    label=leadfield.label;

    save(fullfile(save_dir, 'SPMgainmatrix_cord1.mat'), 'Gx', 'label', spm_get_defaults('mat.format'));

    save(fullfile(save_dir, 'SPMgainmatrix_cord2.mat'), 'Gy', 'label', spm_get_defaults('mat.format'));

    save(fullfile(save_dir, 'SPMgainmatrix_cord3.mat'), 'Gz', 'label', spm_get_defaults('mat.format'));


    % PREP SENSOR CHANNEL INFO
    chanind = D.indchantype('MEG');
    channames = D.chanlabels(chanind);
    sensors = D.sensors('MEG');

    chanpos_s = zeros(length(chanind),3);
    chanori_s = zeros(length(chanind),3);
    headind = []; backind = []; frontind = [];

    for f = 1:length(chanind)
        idx = find(strcmp(sensors.label, channames{f}));
        chanpos_s(f,:) = sensors.coilpos(idx,:);
        chanori_s(f,:) = sensors.coilori(idx,:);

        if chanpos_s(f,2) > 200
            headind = [headind chanind(f)];
        elseif chanpos_s(f,2) < 200 && chanpos_s(f,3) < 0
            backind = [backind chanind(f)];
        elseif chanpos_s(f,2) < 200 && chanpos_s(f,3) > 0
            frontind = [frontind chanind(f)];
        end
    end

    % BEAMFORMER PROCESSING LOOP

    for nn = 1% (runs)
        Dbf = D;
        chanind = indchantype(Dbf, 'MEGMAG');
        chanind = chanind(~ismember(chanind, Dbf.badchannels));
        channames = Dbf.chanlabels(chanind);
        sensors = Dbf.sensors('MEG');

        chanpos_s = zeros(length(chanind), 3);

        for f = 1:length(chanind)
            idx = find(strcmp(sensors.label, channames{f}));
            chanpos_s(f,:) = sensors.coilpos(idx,:);
        end

        bodyind = chanind(chanpos_s(:,2) < 200 & chanpos_s(:,3) < 0);
        covsamples = find(Dbf.time > covtime(1) & Dbf.time < covtime(2));
        Ntrials = Dbf.ntrials;
        allDbf=Dbf(:,:,:);

        % Covariance matrix
        covdata = zeros(length(bodyind));
        for f = 1:Ntrials
            data = squeeze(allDbf(bodyind, covsamples, f));
            data = data - mean(data, 2);
            covdata = covdata + data * data';
        end
        covdata = covdata / Ntrials;
        invcovdata = pinv(covdata + lambda * trace(covdata) * eye(size(covdata)));

        % Load lead fields for body channels
        lfbodyind = find(ismember(leadfield.label, Dbf.chanlabels(bodyind)));
        weights = zeros(Ncordpoints, length(bodyind), 3);

        L=zeros(num_sources,numel(leadfield.label),3);
        L(:,:,1)=Gx;
        L(:,:,2)=Gy;
        L(:,:,3)=Gz;%just do this 3d for now

        for cpind = 1:Ncordpoints
            H = squeeze(L(cpind, lfbodyind,:));
            weights(cpind, :, :) = (pinv(H' * invcovdata * H) * H' * invcovdata)';
        end

        % Apply beamformer
        wdata = zeros(Ncordpoints, Ntrials, length(covsamples), 3);
        tcord = zeros(Ncordpoints, length(covsamples), 3);

        for cpind = 1:Ncordpoints
            for f = 1:Ntrials
                data = squeeze(allDbf(bodyind, covsamples, f));
                data = data - mean(data, 2);
                wdata(cpind, f, :, :) = data' * squeeze(weights(cpind, :, :));
            end

        end
        stidx=find(strcmp(D.conditions,'static'));
        reidx=find(strcmp(D.conditions,'rest'));

        wdatarest=wdata(:,reidx,:,:);
        wdatastatic=wdata(:,stidx,:,:);

        % Save results
        if HFC
            bffilename = fullfile(save_dir, ['bfdata_', sub]);

        else
            bffilename = fullfile(save_dir, ['bfdata_nohfc_', sub]);
        end
        spmfilename = Dbf.fullfile;
        save(bffilename, 'wdata', 'weights', 'covsamples', 'lambda', ...
            'spmfilename', 'bodyind', 'headind', 'backind', 'frontind', 'wdatarest', 'wdatastatic')%,'brainsrc','headmesh','source_brain');


    end
end
%now we have time series for each source point (only good channels)

