%source freq analysis and coherence

close all
clear all
clc

sub='OP00226';
freqband=[10 35];
which_ori='all'; %'all', '1', '2', '3'

addpath('D:\brainspineconnectivity\stats')


%% filenames for brain spinal cord and emg------------------------------------------

save_dir = fullfile('D:\MSST001', [sub '_contrast']);
if ~exist(save_dir, 'dir')
    mkdir(save_dir)
end
bffilename = fullfile(save_dir, ['bfdata_', sub]);


datwithEMGmerged = fullfile('D:\MSST001', ...
    ['sub-' sub], ...
    'ses-001', ...
    'meg', ...
    ['pmergedoe1000mspddfflo45hi45hfcstatic_001_array1.mat']);

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
%% t tests:

n_permutations = 100;

tvals = zeros(1, num_sources);  
pvals = zeros(1, num_sources);


    for k = 1:num_sources
        statdat = log(powstat(1:ntrials, k));
        restdat = log(powrest(1:ntrials, k));

        [~, p, ~, stats] = ttest2(statdat, restdat, 'Tail', 'left');
        warning('one-tailed, testing for power decrease')
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
    %set(f,'Position',[693,1072,849,292])
    set(gca, 'YDir', 'reverse') %neg upwards
    scatter(src.pos(peakind,2),pk,'r*')

%correct p-values
[~, ~, ~, adj_p]=fdr_bh(pvals,0.05,'dep','no');
pSig=adj_p<0.05;

% if adj_p(peakind)
%     disp('Sig peak')
% else
%     disp('NOT SIG')
% end

%get source orientation
%getOptOri

plot_func_brainspine_smooth(mtorso, src, tvals, [], [], 'raw')
% quiver3(src.pos(peakind,1), src.pos(peakind,2), src.pos(peakind,3), ...
%         vecs(peakind,1), vecs(peakind,2), vecs(peakind,3), ...
%         25, 'LineWidth', 2, 'MaxHeadSize',3,'Color', [1 0.4 0]);



savename=sprintf('spine_source_%s_%s',sub,which_ori);
save(fullfile(save_dir,savename), 'pSig', 'tvals','src','thresholds','peakind')

%% extract power from each trial and correlate:

%corr between spinal cord and emg
% sc=log(powstat(1:ntrials,peakind2use));
% emg=log(mean(squeeze(fstat.powspctrm(1:ntrials,end,:)),2)); % tr X freq

%load brain bf data
% braindat = load(fullfile(save_dir,['chandata_brain_', sub]));
% 
% brain=log(braindat.pkData(1:ntrials));

% brainchanidx=find(contains(grad.label,braindat.brainlabmax));
% maxbrainpos=grad.chanpos(brainchanidx,:);




% srcbrain=braindat.brainsrc;
% inside_br = srcbrain.pos(srcbrain.inside, :);
% plot3(inside_br(maxidx,1), inside_br(maxidx,2),inside_br(maxidx,3), ...
%     'o', 'MarkerSize', 10, 'MarkerFaceColor', 	[0.5 0 0.8], 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);


% plot3(src.pos(peakind2use,1), src.pos(peakind2use,2),src.pos(peakind2use,3),'ro','MarkerSize',16)
% plot3(braindat.brainsrc.pos(braindat.idx2plot,1), braindat.brainsrc.pos(braindat.idx2plot,2),braindat.brainsrc.pos(braindat.idx2plot,3),'ro','MarkerSize',16)
% 

%plot_func_brainspine_dat(mtorso,src,tvals,braindat.brainsrc, braindat.func,grad)


% idx2del=find(emg>3);
% brain(idx2del)=[]; sc(idx2del)=[]; emg(idx2del)=[];

%[r, p] = corr(sc, brain, 'Type', 'Spearman');
% [r, p] = corr(emg, brain, 'Type', 'Spearman');
% [r, p] = corr(sc, emg, 'Type', 'Spearman');


% figure; plot(sc,brain,'k.'); lsline
% title('spinal cord and brain')
% box off
% 
% text_x = min(sc) + 0.05 * range(sc);  % X position for text
% text_y = max(brain) - 0.1 * range(brain);  % Y position for text
% 
% % Format the string with r and p
% txt = sprintf('r = %.2f, p = %.2g', r, p);
% 
% % Add the text to the plot
% text(text_x, text_y, txt, 'FontSize', 12, 'FontWeight', 'bold');
% set(gca, 'FontSize', 12);  % Applies to X and Y tick labels
% 
% 
% % figure; plot(emg,brain,'k.'); lsline
% title('emg and brain')
% 
% figure; plot(emg, sc, 'k.'); lsline
% title('emg and sc')


%extract time series:
% cfg=[];
% cfg.channel=[alldat_stat.label(peakind_pos);alldat_stat.label(end)];
% scandemg=ft_selectdata(cfg,alldat_stat);
% 
% cfg=[];
% cfg.channel=braindat.brainlabmax;
% brain=ft_selectdata(cfg,dat_static);
% 
% allts=ft_appenddata([],brain,scandemg);
% 
% 
% 
%     cfg            = [];
%     cfg.output     = 'fourier';
%     cfg.method     = 'mtmfft';
%     cfg.foilim     = [0 40];
%     cfg.tapsmofrq  = '2';
%     cfg.keeptrials = 'yes';
%    % cfg.pad='maxperlen';
% 
%     freqall   = ft_freqanalysis(cfg, allts);


%     cfg=[];
% 
%     cfg.method='powcorr_ortho';
%     cfg.removemean='yes';
%     coh=ft_connectivityanalysis(cfg,freqall);
% 
%     figure;
%     box off
%     subplot(131)
%     plot(coh.freq, abs(squeeze(coh.powcorrspctrm(1,2,:)))); %brain
%     title('brain spine')
%     ylim
%      set(gca, 'FontSize', 12);
%     subplot(132)
%     plot(coh.freq, abs(squeeze(coh.powcorrspctrm(1,3,:))));
%     title('brain emg')
%     subplot(133)
%     plot(coh.freq, abs(squeeze(coh.powcorrspctrm(2,3,:))))
%     title('spine emg')
%     set(gca, 'FontSize', 12);  % Applies to X and Y tick labels


    %%
%     cfg=[];
% 
%     cfg.method='coh';
%     cfg.complex='absimag';
%     coh=ft_connectivityanalysis(cfg,freqall);
% 
%     figure; hold on
%     subplot(131)
%     plot(coh.freq, squeeze(coh.cohspctrm(1,2,:))); %brain
%     title('brain  spine')
%     subplot(132)
%     plot(coh.freq, squeeze(coh.cohspctrm(1,3,:)));
%     title('brain  emg')
%     subplot(133)
%     plot(coh.freq, squeeze(coh.cohspctrm(2,3,:)))
%     title('spine  emg')
% 
% 
% savename=sprintf('cohall_%s',sub);
% save(fullfile(save_dir,savename), 'coh')

%%



% cfg            = [];
% cfg.output     = 'fourier';
% cfg.method     = 'mtmfft';
% cfg.foilim     = [5 75];
% cfg.tapsmofrq  = 5;
% cfg.keeptrials = 'yes';
% 
% fstat   = ft_freqanalysis(cfg,alldat_stat);
% frest   =ft_freqanalysis(cfg,alldat_rest);
% 
% seed='EXG1';
% spine=alldat_stat.label(1:end-1);
% tmpcmb=repmat({seed}, length(spine),1);
% combs=[spine tmpcmb];
% 
% cfg            = [];
% cfg.method     = 'coh';
% %cfg.complex    ='absimag';
% cfg.channelcmb = combs;
% coh      = ft_connectivityanalysis(cfg,fstat);
% cohrest  = ft_connectivityanalysis(cfg,frest);
% 
% coh2plot=mean(coh.cohspctrm,2);
% 
% figure;
% plot_func_spine_dat(mtorso,src,coh2plot,grad)
% plot3(src.pos(peakind,1), src.pos(peakind,2),src.pos(peakind,3),'ro','MarkerSize',12)
% title('coh diff')

%% coherence difference

% N = 549;
% C1=coh.cohspctrm(:,betafreq);
% C2=cohrest.cohspctrm(:,betafreq);
% 
% var1 = (2 .* (1 - C1.^2).^2) ./ N;
% var2 = (2 .* (1 - C2.^2).^2) ./ N;
% 
% % Step 2: Compute Z-score for each matrix element
% Z = (C1 - C2) ./ sqrt(var1 + var2);
% 
% meanZ=mean(Z,2);
% [~,peakind]=max(meanZ); %maxk?
% 
% 
% % cohdiff=coh.cohspctrm(:,betafreq)-cohrest.cohspctrm(:,betafreq);
% % plot_func_spine_dat(mtorso,src,mean(cohdiff,2),grad)
% 
% figure;
% plot_func_spine_dat(mtorso,src,meanZ,grad)
% plot3(src.pos(peakind,1), src.pos(peakind,2),src.pos(peakind,3),'ro','MarkerSize',12)
% title('coh diff')
% 
% figure; plot(coh.freq,coh.cohspctrm(peakind,:))


%CVA (this is spine-EMG)
% betafreq=find(fstat.freq >=15 & fstat.freq <=30);
% Xdat=squeeze(fstat.fourierspctrm(:,peakind2use,betafreq));
% 
% Ydat=squeeze(fstat.fourierspctrm(:,end,betafreq));
% 
% usefreq=fstat.freq (betafreq);
%CVA_coh(Ydat,Xdat,usefreq,'emgspine',sub)


% X2=mean(fstat.powspctrm(1:end-1,:),2)';
% Y2=mean(fstat.powspctrm(end,:),2)';
% d=[33 1];
% [V,rho]=mcca([X2 Y2],[33 1])

% Find where pSig == 1
% isOne = pSig == 1;
% 
% % Detect transitions
% d = diff([0 isOne 0]);  % Pad to catch edges
% 
% starts = find(d == 1);
% ends   = find(d == -1) - 1;

% Extract valid clusters
% clusters = {};
% cluster_sums = [];  % Store sum of abs t-values for each cluster
% 
% for i = 1:length(starts)
%     cluster = starts(i):ends(i);
%     if numel(cluster) > 1
%         t_cluster = tvals(cluster);
%         % Check if all t-values are of the same sign
%         if all(t_cluster > 0) || all(t_cluster < 0)
%             clusters{end+1} = cluster;
%             cluster_sums(end+1) = sum(abs(t_cluster));
%         end
%     end
% end

% Find cluster with maximum summed absolute t-stat
% if ~isempty(cluster_sums)
%     [~, maxIdx] = max(cluster_sums);
%     maxCluster = clusters{maxIdx};
% else
%     maxCluster = [];  % No valid clusters found
%     fprintf('No valid sig clusters\n')
% end


% figure;
% plot_func_spine_dat(mtorso,src,tvals,grad)
% if ~isempty(maxCluster)
%     if all(tvals(maxCluster) > 0)
%         fprintf('Max cluster has positive t-values.\n');
%     elseif all(tvals(maxCluster) < 0)
%         fprintf('Max cluster has negative t-values.\n');
%     end
%     plot3(src.pos(maxCluster,1), src.pos(maxCluster,2),src.pos(maxCluster,3),'ro','MarkerSize',16)
% end
% 
% title('max cluster')

