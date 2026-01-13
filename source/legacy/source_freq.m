%source freq analysis

close all
clear all
clc

subs = {'OP00212', 'OP00213', 'OP00215', 'OP00219', ...
    'OP00220', 'OP00221', 'OP00225', 'OP00226'};

%subs={'OP00224'}; %_002


freqband=[10 35];
ori2test={'all'};%, '1', '2', '3'};
HFC=1;
addpath('D:\brainspineconnectivity\stats')


%% filenames for brain spinal cord and emg------------------------------------------
for s=1:length(subs)

    sub=subs{s};

    for oo=1:length(ori2test)
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
                ['pmergedoe1000mspddfflo45hi45hfcstatic_001_array1.mat']);
        else
            bffilename = fullfile(save_dir, ['bfdata_nohfc_', sub]);
            datwithEMGmerged = fullfile('D:\MSST001', ...
                ['sub-' sub], ...
                'ses-001', ...
                'meg', ...
                ['pmergedoe1000mspddfflo45hi45static_001_array1.mat']);
        end


        geomfile = 'D:\MSST001\generic_merged\geoms.mat';


        %% load
        D=spm_eeg_load(datwithEMGmerged);
        grad=D.sensors('MEG');

        castforward = load(geomfile);
        mtorso = castforward.mesh_torso;
        src=castforward.sources_center_line;

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

        powstat=zeros(D.ntrials,num_sources); %save power for each source pt

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
        %% t tests and permutation testing

        n_permutations = 100;

        tvals = zeros(1, num_sources);
        pvals = zeros(1, num_sources);


        for k = 1:num_sources
            statdat = log(powstat(1:ntrials, k));
            restdat = log(powrest(1:ntrials, k));

            [~, p, ~, stats] = ttest2(statdat, restdat, 'Tail', 'left');
            warning('one-tailed dep test, testing for power decrease') %i.e. neg t vals
            tvals(k) = stats.tstat;
            pvals(k) = p;
        end

        null_tvals = zeros(num_sources, n_permutations);

        for i = 1:n_permutations
            for k = 1:num_sources
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

 
        f=figure; plot(src.pos(:,2),tvals, 'ko-');
        hold on; plot(src.pos(:,2), thresholds,'k--','LineWidth',1)
        box off
        set(gca,'FontSize',14)
        set(f,'Position',[693,1072,849,292])
        set(gca, 'YDir', 'reverse') %neg upwards
        scatter(src.pos(peakind,2),pk,'r*')
        title(sprintf('subject %s',sub))

        %correct p-values
%         [~, ~, ~, adj_p]=fdr_bh(pvals,0.05,'dep','no');
%         pSig=adj_p<0.05;



        %get source orientation
        %getOptOri

        %plot_func_brainspine_smooth(mtorso, src, tvals, [], [], 'raw')
        % quiver3(src.pos(peakind,1), src.pos(peakind,2), src.pos(peakind,3), ...
        %         vecs(peakind,1), vecs(peakind,2), vecs(peakind,3), ...
        %         25, 'LineWidth', 2, 'MaxHeadSize',3,'Color', [1 0.4 0]);


            static = log(powstat(1:ntrials, :));
            rest = log(powrest(1:ntrials, :));


%         if HFC
%             savename=sprintf('spine_source_%s_%s',sub,which_ori);
%         else
%             savename=sprintf('spine_source_%s_%s_nohfc',sub,which_ori);
% 
%         end
%         save(fullfile(save_dir,savename), 'pSig', 'tvals','src','thresholds','peakind', 'static', 'rest')
     end
 end
