%source freq analysis for classical group stats

close all
clear all
clc

%subs = {'OP00212', 'OP00213', 'OP00215', 'OP00219', ...
%'OP00220', 'OP00221', 'OP00225', 'OP00226'};

subs={'OP00224'}; %_002

HFC=1;
freqband=[10 35];
ori2test={'all', '1', '2', '3'};

geomfile = 'D:\MSST001\generic_merged\geoms.mat';
castforward = load(geomfile);
src=castforward.sources_center_line;


for s=1:length(subs)
    for oo=1:length(ori2test)

        sub=subs{s};
        which_ori=ori2test{oo};
        save_dir = fullfile('D:\MSST001', [sub '_contrast']);

        if ~exist(save_dir, 'dir')
            mkdir(save_dir)
        end

        if HFC
            bffilename = fullfile(save_dir, ['bfdata_', sub]);


            datwithEMGmerged = fullfile('D:\MSST001', ...
                ['sub-' sub], ...
                'ses-001', ...
                'meg', ...
                ['pmergedoe1000mspddfflo45hi45hfcstatic_002_array1.mat']);

        else
            bffilename = fullfile(save_dir, ['bfdata_nohfc_', sub]);

            datwithEMGmerged = fullfile('D:\MSST001', ...
                ['sub-' sub], ...
                'ses-001', ...
                'meg', ...
                ['pmergedoe1000mspddfflo45hi45static_001_array1.mat']);


        end


        %% load
        D=spm_eeg_load(datwithEMGmerged);
        grad=D.sensors('MEG');

        st=load(bffilename);

        dat=spm2fieldtrip(D);

        %crop struct to n channels = n sources
        num_sources=size(st.wdatastatic,1);

        %structure to put source data into
        cfg=[];
        cfg.channel = dat.label(1:num_sources);
        alldat=ft_selectdata(cfg,dat);

        %process EMG
        cfg=[];
        cfg.channel='EXG1';
        cfg.rectify ='yes';
        cfg.detrend='yes';
        EMG=ft_preprocessing(cfg, dat);

        alldat_EMG=ft_appenddata([], alldat, EMG);


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

            for k = 1:numel(alldat.trial)
                alldat.trial{k}(:,:) = squeeze(st.wdata(:,k,:,o)); %alldat does not have EMG
            end


            if o == ori_list(1)
                for j = 1:numel(alldat.label)-1
                    alldat.label{j} = sprintf('source%g', j);
                end
            end

            % Spectral power analysis
            cfg = [];
            cfg.output     = 'pow';
            cfg.method     = 'mtmfft';
            cfg.foilim     = freqband;
            cfg.tapsmofrq  = 2;
            cfg.keeptrials = 'yes';
            cfg.pad        = 'nextpow2';

            fdat = ft_freqanalysis(cfg, alldat);

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
        %% calculate contrast

        % meanstattr=mean(powstat,1); %mean over trials
        % meanresttr=mean(powrest,1);

        ntrials=min([size(powstat,1) size(powrest,1)]);
        powdiff=log(powrest(1:ntrials,:))-log(powstat(1:ntrials,:)); %diff for each source point and trial; reduction produces positive values


nBins = 50;
data=powdiff;
% Define bin edges across all data
allDataMin = min(data(:));
allDataMax = max(data(:));
binEdges = linspace(allDataMin, allDataMax, nBins + 1);

% Initialize histogram matrix (bins x sourcepoints)
histMatrix = zeros(nBins, size(data, 2));

% Loop through each sourcepoint and compute histogram
for i = 1:size(data, 2)
    histCounts = histcounts(data(:, i), binEdges);
    histMatrix(:, i) = histCounts';
end

% Plot heatmap
figure;
imagesc(histMatrix);
colormap('hot');       % or 'parula', 'viridis' (if you have it), etc.
colorbar;
xlabel('Sourcepoint');
ylabel('Histogram Bin');
title('Histogram Heatmap of Each Sourcepoint');

% Optional: improve Y-tick labels to reflect bin centers
binCenters = (binEdges(1:end-1) + binEdges(2:end)) / 2;
yticks = round(linspace(1, nBins, 5));  % e.g. 5 tick labels
yticklabels = arrayfun(@(x) sprintf('%.2f', binCenters(x)), yticks, 'UniformOutput', false);
set(gca, 'YTick', yticks, 'YTickLabel', yticklabels);
[~, zeroIdx] = min(abs(binCenters - 0));
hold on;
yline(zeroIdx, 'c-', 'LineWidth', 1.5);

        %[pk, peakind] = max(powdiff);

        %     f=figure; plot(src.pos(:,2),powdiff, 'ko-');hold on
        %     box off
        %     set(gca,'FontSize',14)
        %     %set(f,'Position',[693,1072,849,292])
        %     %set(gca, 'YDir', 'reverse') %neg upwards
        %     scatter(src.pos(peakind,2),pk,'r*')


        if HFC
            savename=sprintf('spine_source_group_%s_%s_fixed',sub,which_ori);
        else
            savename=sprintf('spine_source_group_%s_%s_nohfc_fixed',sub,which_ori);

        end
        save(fullfile(save_dir,savename), 'powdiff')
    end
end









