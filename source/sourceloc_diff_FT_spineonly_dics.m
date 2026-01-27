
% SPINE EMG COH

clear all
close all
clc

addpath('C:\Users\mspedden\Documents\brainspineconnectivity\source')
addpath('C:\Users\mspedden\Documents\spm')
spm('defaults','EEG')

addpath('C:\Users\mspedden\Documents\fieldtrip')
ft_defaults

save_dir='C:\Users\mspedden\Documents\brainspine_save';
rng(1) %for permutation testing

%n=9 for spinal cord analyses
subs = {'OP00212','OP00213',  'OP00215', 'OP00219', ...
    'OP00220', 'OP00221', 'OP00224', 'OP00225', 'OP00226'};

generic_dir = 'C:\Users\mspedden\Documents\new_leadfields_and_geom'; %where I have saved folder with brainspine leadfields and geoms (meshes)
geomfile = fullfile(generic_dir, 'geometries_cervical_realistic.mat');

LFop='spine'; %only want leadfields from spine here.
rectify=1; %EMG
fband=[10 35];
mult_comp_corr=1;

subjResults=struct();

for ss=1:length(subs)

    sub=subs{ss};

    if strcmp(sub, 'OP00224') %saved under 002
        datwithEMGmerged = fullfile('C:\Users\mspedden\Documents', ...
            ['sub-' sub], ...
            'ses-001', ...
            'meg', ...
            'pmergedoe1000mspddfflo45hi45hfcstatic_002_array1.mat');

    else

        datwithEMGmerged = fullfile('C:\Users\mspedden\Documents', ...
            ['sub-' sub], ...
            'ses-001', ...
            'meg', ...
            'pmergedoe1000mspddfflo45hi45hfcstatic_001_array1.mat');
    end

    load(geomfile)
    D=spm_eeg_load(datwithEMGmerged);
    grad_mm=D.sensors('MEG');
    ftdat = spm2fieldtrip(D);

    badchans=D.chanlabels(D.badchannels);

    %remove bad channels here.
    cfg=[];
    cfg.channel=setdiff(ftdat.label,badchans);
    ftdat=ft_selectdata(cfg,ftdat);

    %% rectify EMG
    if rectify
        cfg=[];
        cfg.rectify='yes';
        cfg.channel='EXG1';
        ftdatr=ft_preprocessing(cfg,ftdat);

        for k=1:length(ftdat.trial)
            ftdat.trial{k}(end,:)=ftdatr.trial{k}; %ftdat has rectified emg
        end

    end
%% load and organise spinal cord leadfields
    [Gx, Gy, Gz] = build_leadfield_matrices(fullfile(generic_dir,'cervical_realistic_brain_spine'), LFop);

    nsourcepoints = size(Gx,1);
    nchannels     = size(Gx,2);
    spchanidx=find(grad_mm.coilpos(:,2) < 200); %indexed locally in grad struct (same indexing as LF)
    spchanlabs=grad_mm.label(spchanidx);

    %% clip leadfields to spinal cord channels only
    Gx=Gx(:,spchanidx);
    Gy=Gy(:,spchanidx);
    Gz=Gz(:,spchanidx);
    
    %put leadfields into fieldtrip format
    Lf.pos    = sources_cent.pos;     % nsourcepoints x 3
    Lf.inside = sources_cent.inside;    
    Lf.unit   = 'mm';
    Lf.label  = grad_mm.label(spchanidx);   % nchannels x 1 cell
    Lf.leadfielddimord = '{pos}_chan_ori';
    Lf.leadfield = cell(1,nsourcepoints);

    for k = 1:nsourcepoints
        % Combine X/Y/Z components like FT is used to
        Lf.leadfield{k} = [Gx(k,:)' Gy(k,:)' Gz(k,:)']; % nchannels x 3
    end

    % 2. dummy head model for input config only (not actually used)
    cfg                     = [];
    cfg.method              = 'infinite';
    cfg.siunits=1;
    cfg.grad=grad_mm;
    cfg.conductivity = 1;

    dummyvol = ft_prepare_headmodel(cfg,mesh_torso);


    %% beamforming----------------------------
    %1. get trial wise freq dat
    cfg=[];
    cfg.output     = 'powandcsd';
    cfg.method     = 'mtmfft';
    cfg.foilim     = fband;
    cfg.tapsmofrq  = 1;
    cfg.keeptrials='yes';
    freqdat_tr=ft_freqanalysis(cfg,ftdat); 

    cfg=[];
    cfg.avgoverfreq='yes';
    freqdat_tr=ft_selectdata(cfg,freqdat_tr);%need to average over frequency - required input to permutation test

    %% separate conditions

    statidx=find(ftdat.trialinfo==1);
    restidx=find(ftdat.trialinfo==2);

    [nTrials,~] = min([length(statidx) length(restidx)]); %ensure same n trials across conditions

    cfg=[];
    cfg.trials=statidx(1:nTrials);
    statdat=ft_selectdata(cfg,freqdat_tr); %trial wise data pr condition

    cfg.trials=restidx(1:nTrials);
    restdat=ft_selectdata(cfg,freqdat_tr);

    %% get struct with mean freq data for constructing bf filter
    cfg=[];
    cfg.output     = 'powandcsd';
    cfg.method     = 'mtmfft';
    cfg.foilim     = fband;
    cfg.tapsmofrq  = 1;
    cfg.keeptrials = 'no';
    freqdat=ft_freqanalysis(cfg,ftdat);

    cfg=[];
    cfg.avgoverfreq='yes';
    freqdat=ft_selectdata(cfg,freqdat);

    %% DICS filter constructed on both conditions

    cfg=[];
    cfg.grid = sources_cent;
    cfg.headmodel=dummyvol;
    cfg.sourcemodel.leadfield=Lf;
    cfg.dics.keepfilter='yes';
    cfg.dics.lambda=10;
    cfg.method = 'dics';
    cfg.refchan='EXG1'; %EMG

    coh_source = ft_sourceanalysis(cfg,freqdat);

    %% apply to each condition separately, do permutation test
    cfg=[];
    cfg.grid = sources_cent;
    cfg.headmodel=dummyvol;
    cfg.sourcemodel.leadfield=Lf;
    cfg.dics.filter=coh_source.avg.filter;
    cfg.dics.lambda=10;
    cfg.method = 'dics';
    cfg.refchan='EXG1';

    source_stat = ft_sourceanalysis(cfg,statdat);
    source_rest = ft_sourceanalysis(cfg,restdat);

    cfg.permutation = 'yes';
    cfg.numpermutation=500;
    source_perm = ft_sourceanalysis(cfg, statdat, restdat); %permutation test

    nPerm = numel(source_perm.trialA); 
    cohDiff_perm = zeros(nsourcepoints, nPerm);

    for i = 1:nPerm
        cohA_perm = source_perm.trialA(i).coh;
        cohB_perm = source_perm.trialB(i).coh;
        cohDiff_perm(:, i) = cohA_perm - cohB_perm; %only interested in contraction > rest
    end

    maxPerm = max(cohDiff_perm, [], 1); % max over sources

    if mult_comp_corr
        thr95 = prctile(maxPerm, 95,2); %1 x 1 global thresh
    else
        thr95=prctile(cohDiff_perm,95,2); %95th percentile for each source over all permutations
        warning('Uncorrected thresh might not work with plots')
    end

%% source structure with coh difference

    coh_diff=source_perm.avgA.coh-source_perm.avgB.coh;
    source_diff=coh_source;
    source_diff.avg.coh=coh_diff;

     % One-sided permutation p-value 
pvals = zeros(nsourcepoints,1);
for s = 1:nsourcepoints
    permDist = cohDiff_perm(s, :);
    obsVal   = coh_diff(s);
    pvals(s) = (sum(permDist >= obsVal) + 1) / (nPerm + 1);
end
invp = -log10(pvals);

mask = coh_diff > thr95;    % same as before

% Mask the inverse-p map
invp_masked = invp;
invp_masked(~mask) = 0;     % or NaN if


    % sources structure with p values
    source_p=coh_source;
    source_p.avg.coh=invp_masked;

    cfg = [];
    cfg.parameter = 'coh';
    spine_int = ft_sourceinterpolate(cfg, source_p, mesh_wm);

%     permMean = mean(cohDiff_perm, 2);
%     permStd  = std(cohDiff_perm, 0, 2);
%     
%     z_coh = (coh_diff - permMean) ./ permStd;
%     permStd(permStd == 0) = NaN;
% 
%     source_z = coh_source;        %source structure
%     source_z.avg.coh = z_coh; 

    %interpolate on spinal cord mesh
%     cfg = [];
%     cfg.parameter = {'coh'};
%     cfg.interpmethod='nearest';
%     spine_int=ft_sourceinterpolate(cfg,source_diff, mesh_wm);
% 
%     cfg = []; %interpolate z score
%     cfg.parameter = {'coh'};
%     spine_int_z=ft_sourceinterpolate(cfg,source_z, mesh_wm);
%     

    %% add a significance mask
    source_mask=coh_source; %copy
    source_mask.avg.pow = coh_diff > thr95;
    maxidx=NaN;

    %before interpolating mask find sig point where diff is max
    %this is all for looking at orientation of source where sig diff is max
    sigMask=source_mask.avg.pow;
    if any(sigMask) %if any sig points
        maskedDiff = coh_diff;
        maskedDiff(~sigMask) = -inf;   % remove nonsignificant
        [~, maxidx] = max(maskedDiff); %max idx is nan if no sig diffs
    end
% 
%     fprintf('Selected source point index: %d\n', maxidx);


source_mask = coh_source;
source_mask.avg.coh = double(mask);  
cfg.parameter = 'coh';
cfg.interpmethod = 'nearest';
source_mask_int = ft_sourceinterpolate(cfg, source_mask, mesh_wm);
spine_int.mask = source_mask_int.coh;  % same as before

% 
%     if isnan(maxidx)
%         subjResults(ss).maxdiff.pos = [NaN NaN NaN];
%         subjResults(ss).maxdiff.ori = [NaN NaN NaN];
%     end
% 
% Lf_voxel = Lf.leadfield{maxidx};      
% filter_voxel = coh_source.avg.filter{maxidx};  
% pos=Lf.pos(maxidx,:);
% chanIdx = ismember(freqdat.label, coh_source.avg.label); 
% 
% %% get csd
% nCh=sum(chanIdx);
% C = zeros(nCh, nCh);
% 
% % build the full matrix manually from labelcmb and crsspctrm
% for k = 1:length(freqdat.crsspctrm)
%     ch1 = find(strcmp(freqdat.label, freqdat.labelcmb{k,1}));
%     ch2 = find(strcmp(freqdat.label, freqdat.labelcmb{k,2}));
%     if chanIdx(ch1) && chanIdx(ch2)
%         idx1 = find(find(chanIdx) == ch1);
%         idx2 = find(find(chanIdx) == ch2);
%         C(idx1, idx2) = freqdat.crsspctrm(k);
%         C(idx2, idx1) = conj(freqdat.crsspctrm(k)); % Hermitian
%     end
% end
% 
% P = filter_voxel * C * filter_voxel.';   % 3 x 3
% 
% % Eigendecomposition to get dominant orientation
% [V, D] = eig(P);                        % eigenvectors & eigenvalues
% [~, idx] = max(diag(D));                % index of largest eigenvalue
% 
% %% Dominant dipole orientation (normalized)
% dip_orient = V(:, idx);
% dip_norm = dip_orient / norm(dip_orient);

% figure; hold on;
%     scatter3(sources_cent.pos(:,1), sources_cent.pos(:,2), sources_cent.pos(:,3), 5, 'k');
%     quiver3(pos(1), pos(2), pos(3), dip_norm(1), dip_norm(2), dip_norm(3), 8, 'r', 'LineWidth', 2);
%     xlabel('x'); ylabel('y'); zlabel('z'); axis equal;
%     title(sprintf('Dipole Orientation at Voxel %d (Contraction Condition)', maxidx));
%     grid on;

%     %interpolate the mask
%     cfg = [];
%     cfg.parameter = 'pow';
%     cfg.interpmethod = 'nearest';
%     source_mask_int = ft_sourceinterpolate(cfg, source_mask, mesh_wm);
%     spine_int.mask=source_mask_int.pow; %add to struct for plotting
%     spine_int_z.mask=source_mask_int.pow;

    %clip torso mesh to plot less of lower body, should save this
    y=mesh_torso.vertices(:,2);
    keep_vert=y>-61;
    new_idx=zeros(size(keep_vert));
    new_idx(keep_vert)=1:sum(keep_vert);
    faces_keep=all(keep_vert(mesh_torso.faces),2);
    new_faces=new_idx(mesh_torso.faces(faces_keep,:));
    new_vertices=mesh_torso.vertices(keep_vert,:);
    mesh_cut.vertices=new_vertices;
    mesh_cut.faces=new_faces;
    mesh_cut.unit=mesh_torso.unit;

    % source plots

    ncol = 256;
    addpath('C:\Users\mspedden\Documents\fieldtrip\external\matplotlib\')
    brain_color = [0.92 0.92 0.92];   % first color = background
    hotmap = flipud(magma(ncol-1));
    cmap = [brain_color; hotmap];

    figure
    cfg = [];
    cfg.figure='gcf';
    cfg.method = 'surface';
    cfg.funparameter = 'coh';
    cfg.funcolormap =cmap;
    cfg.funcolorlim=[0 max(spine_int.coh(:))];
    cfg.projmethod = 'nearest';
    cfg.surffile = mesh_wm;
    ft_sourceplot(cfg, spine_int);
    view( -250, -1)
    camlight
    ax = gca;                 
    ax.FontSize = 14;
   hpatch = findobj(gcf, 'Type', 'patch');
    set(hpatch, 'FaceAlpha',0.9)

%srcmax = sources_cent.pos(maxidx,:);

% % compute nearest vertex on the surface mesh
% d = vecnorm(mesh_wm.vertices - srcmax, 2, 2);
% [~, nearestIdx] = min(d);
% xyz_surf = mesh_wm.vertices(nearestIdx, :);
% 
% hold on;
% plot3(xyz_surf(1), xyz_surf(2), xyz_surf(3), 'o', ...
%       'MarkerSize', 12, 'MarkerFaceColor', 'cyan');

%dir = [0 0 1];       % arrow pointing upward in Z
%dir = dir ./ norm(dir);

%arrowLength = 15;  % length of arrow

% Base of arrow = head - length * direction
%base = srcmax - dir * arrowLength;

% Plot arrow
% quiver3( base(1), base(2), base(3), ...   
%          dir(1)*arrowLength, dir(2)*arrowLength, dir(3)*arrowLength, ... % vector
%          0, ...  % do not auto-scale
%          'LineWidth', 2, ...
%          'Color', 'black', ...
%          'MaxHeadSize', 30 );

%% save results

    subjResults(ss).coh_diff=coh_diff;        % source_perm.avgA.coh - avgB.coh
    subjResults(ss).thr95=thr95;
    subjResults(ss).sig_mask = spine_int.coh; 
    subjResults(ss).pos = spine_int.pos; 
    subjResults(ss).inside = spine_int.inside;
    subjResults(ss).rest =source_rest; %prob dont need this
    subjResults(ss).stat=source_stat;
    %subjResults(ss).maxdiff.pos=pos; %dont think I need this
    %subjResults(ss).maxdiff.ori=dip_norm;
end
save('groupRes_spine_DICS.mat', 'subjResults')
nSubjects = length(subjResults);
sig_pos = false(nSubjects,1);

%number of subjects with at least one significant source anywhere
for ss = 1:nSubjects
    diffCoh = subjResults(ss).coh_diff;   % source_perm.avgA.coh - avgB.coh
    thr95    = subjResults(ss).thr95; % 95th percentile from permutation

    if any(diffCoh > thr95)
        sig_pos(ss) = true;
    end
end

fprintf('Permutation: %d/%d subjects show a positive effect above threshold\n', ...
    sum(sig_pos), nSubjects);


%% group prevalence

nSubs = length(subjResults);

all_masks = cat(2, subjResults(:).sig_mask);

% Compute prevalence 
group_prevalence = mean(all_masks, 2);

% Create a group source structure for plotting
group_source = subjResults(1); % use one subject as template
group_source.pow = group_prevalence;
group_source = rmfield(group_source, {'coh_diff', 'sig_mask'}); % clean up

group_ft = [];
group_ft.pos = group_source.pos;
group_ft.inside = group_source.inside;
group_ft.pow = group_prevalence;

%% Interpolate group map onto the mesh
threshold = 0.2; %to avoid 'false' overlap interpolation
group_ft.pow(group_ft.pow < threshold) = 0;  % threshold source points

cfg = [];
cfg.parameter = 'pow';
cfg.interpmethod = 'nearest';
group_int = ft_sourceinterpolate(cfg, group_ft, mesh_wm);


% group_int.pow_thresh = group_int.pow;          % copy original
% group_int.pow_thresh(group_int.pow < threshold) = NaN;   % subthreshold = NaN

%% Plot group prevalence map
figure;
cfg = [];
cfg.method = 'surface';
cfg.funparameter = 'pow';
cfg.maskparameter = 'mask';          
cfg.funcolorlim =  [threshold max(group_int.pow)];
cfg.funcolormap = cmap;
cfg.projmethod = 'nearest';
cfg.surffile = mesh_wm;
cfg.opacitylim    = [threshold max(group_int.pow)];
cfg.opacitymap    = 'rampup';
ft_sourceplot(cfg, group_int);
view(90,18);
camlight;
ax = gca;                  
ax.FontSize = 14;
hpatch = findobj(gcf, 'Type', 'patch');
set(hpatch, 'FaceAlpha',0.9)

%% Find contiguous cluster for virtual electrode
all_masks = zeros(nsourcepoints, nSubjects);

for s = 1:nSubjects
    if isempty(subjResults(s).coh_diff)
        mask = zeros(nsourcepoints,1);   % no significant sources
    else
        mask = subjResults(s).coh_diff > subjResults(s).thr95;
        mask = mask(:);             
    end
    all_masks(:,s) = mask;
end

prevalence_loc = mean(all_masks, 2);  % nSourcePoints x 1

mask_thresh = prevalence_loc >= threshold;

% Extract positions of thresholded sources
pos_thresh = sources_cent.pos(mask_thresh,:);

% Compute pairwise distances
distMat = squareform(pdist(pos_thresh));

% Define adjacency points within 6 mm as neighbours
neighborRadius = 6;  
adjMat = distMat < neighborRadius;

G = graph(adjMat);
bins = conncomp(G);   

% Compute cluster sizes
clusterSizes = histcounts(bins, 1:(max(bins)+1));

% Get largest contiguous cluster
[~, idxMax] = max(clusterSizes);

% Indices of points in largest cluster (in pos_thresh)
cluster_idx = find(bins == idxMax);

% Positions of largest cluster in source space
ROIpos = pos_thresh(cluster_idx,:);

%Plot positions of largest cluster
figure; plot(sources_cent.pos(:,2), prevalence_loc)
hold on
for k=1:length(ROIpos)
    plot(ROIpos(k,2), 0.2, 'r*')
end

save(fullfile(save_dir, 'cluster_spineEMG_pos.mat'), 'ROIpos')

%% binomial p and CI
x = sum(sig_pos);
n = nSubjects;
alpha=0.05;
phat = x/n;
lower = betainv(alpha/2, x,     n-x+1);
upper = betainv(1-alpha/2, x+1, n-x);

ci = [lower upper];
p = 1 - binocdf(x-1, n, alpha);



%%
figure;
scatter3(group_ft.pos(:,1), group_ft.pos(:,2), group_ft.pos(:,3), ...
         50, group_prevalence, 'filled');
colorbar;
caxis([0 1]);  % prevalence between 0 and 1
axis equal;
xlabel('X'); ylabel('Y'); zlabel('Z');
title('Group Prevalence (no interpolation)');

%% grave
% figure;

%ft_plot_mesh(mesh_brain, 'facecolor', [0.8 0.3 0.3], 'facealpha', 0.5, 'edgecolor', 'none');
% hold on;
% 
% ft_plot_mesh(mesh_torso, 'facecolor', [0.3 0.3 0.9], 'facealpha', 0.1, 'edgecolor', 'none');
% ft_plot_mesh(mesh_lungs, 'facecolor', [0.8 0.3 0.3], 'facealpha', 1, 'edgecolor', 'none');
% ft_plot_mesh(mesh_heart, 'facecolor', [0.8 0.3 0.3], 'facealpha', 1, 'edgecolor', 'none');
% ft_plot_mesh(mesh_bone, 'facecolor', [0.9 0.85 0.7], 'facealpha', 0.3, 'edgecolor', 'none');
% ft_plot_mesh(mesh_wm, 'facecolor', [0.9 0.9 0.9], 'facealpha', 1, 'edgecolor', 'none');
% ft_plot_sens(ftdat.grad,'coilshape','point','coilsize',8)
% 
% axis equal;
% camlight; lighting gouraud


% figure
% ft_plot_mesh(mesh_torso, 'facecolor', [0.3 0.3 0.9], 'facealpha', 0.1, 'edgecolor', 'none');
% hold on
% plot3(sources_combined.pos(:,1), sources_combined.pos(:,2), sources_combined.pos(:,3),'b.')
%



%% visualise across subjects - 2d

nSubj=numel(subjResults); x=sources_cent.pos(:,2); figure; hold on;

cmap = [
    27,158,119;
    217,95,2;
    117,112,179;
    231,41,138;
    102,166,30;
    230,171,2;
    166,118,29;
    102,102,102;
    55,126,184
    ] / 255;

for s=1:nSubj
    cdiff=subjResults(s).coh_diff; thr=subjResults(s).thr95; sig=cdiff>thr;
    if any(sig), c=cmap(s,:); else, c=[0.7 0.7 0.7]; end
    for i=1:length(x)-1
        if sig(i)&&sig(i+1)
            plot(x(i:i+1),cdiff(i:i+1),'-','Color',c,'LineWidth',1.5,'HandleVisibility','off')
        else
            plot(x(i:i+1),cdiff(i:i+1),'-','Color',[0.7 0.7 0.7],'HandleVisibility','off')
        end
    end
    plot(x(sig),cdiff(sig),'.','Color',c,'MarkerSize',12,'HandleVisibility','off')
    h(s)=plot(nan,nan,'-','Color',c,'LineWidth',1.5);
end
yline(0,':k','HandleVisibility','off')
xlabel('Cranial caudal position (mm)'); ylabel('Coherence difference'); title('Significant coherence differences')
legend(h, arrayfun(@(s) sprintf('Participant %d', s), 1:nSubj, 'UniformOutput', false), ...
    'Location', 'bestoutside');

set(gca, 'FontSize', 13)
grid on;

%% plot orientation max coherence across subjects

% figure; hold on; grid on;
% 
% % define colors for subjects (repeat as needed)
% colors = lines(numel(subjResults));   % nice distinct colors
% 
% scale = 30;  % arrow length scaling factor (adjust to your liking)
% 
% for ss = 1:numel(subjResults)
% 
%     pos = subjResults(ss).maxdiff.pos;    % 1x3
%     ori = subjResults(ss).maxdiff.ori;    % 1x3 or 3x1
% 
%     % ensure row vector
%     ori = ori(:)';  
%     ft_plot_mesh(mesh_wm, ...
%     'facealpha', 0.15, ...            % transparency
%     'facecolor', [0.3 0.6 1.0], ...   % light blue (RGB)
%     'edgecolor', 'none');        hold on
%     quiver3( pos(1), pos(2), pos(3), ...
%              ori(1), ori(2), ori(3), ...
%              scale, ...
%              'Color', colors(ss,:), ...
%              'LineWidth', 2 );
% 
%     % label each subject’s arrow if desired
%     text(pos(1), pos(2), pos(3), sprintf('S%d', ss), ...
%          'Color', colors(ss,:), 'FontSize', 10);
% end
% hold on;  ft_plot_mesh(mesh_brain, ...
%     'facealpha', 0.15, ...            % transparency
%     'facecolor', [0.3 0.6 1.0], ...   % light blue (RGB)
%     'edgecolor', 'none'); 
% 
% xlabel('x'); ylabel('y (up down)'); zlabel('z');
% axis equal;
% title('Ori at max diff source across subjects');

%% sorted by height

% heighttable=readtable('C:\Users\mspedden\Documents\SC_subs_heights.csv');
% heights=heighttable.Var2;
%
% [sortedHeights, sortIdx] = sort(heights, 'descend');  % tallest first
% subjResultsSorted = subjResults(sortIdx);
%
% cmapSorted = cmap(sortIdx, :);
%
% nSubj = numel(subjResultsSorted);
% figure; hold on;
%
% for s = 1:nSubj
%     cdiff = subjResultsSorted(s).coh_diff;
%     thr   = subjResultsSorted(s).thr95;
%     sig   = cdiff > thr;
%
%     if any(sig)
%         c = cmapSorted(s,:);
%     else
%         c = [0.7 0.7 0.7];
%     end
%
%     for i = 1:length(x)-1
%         if sig(i) && sig(i+1)
%             plot(x(i:i+1), cdiff(i:i+1), '-', 'Color', c, 'LineWidth', 1.5, 'HandleVisibility', 'off')
%         else
%             plot(x(i:i+1), cdiff(i:i+1), '-', 'Color', [0.7 0.7 0.7], 'HandleVisibility', 'off')
%         end
%     end
%
%     plot(x(sig), cdiff(sig), '.', 'Color', c, 'MarkerSize', 12, 'HandleVisibility', 'off')
%     h(s) = plot(nan, nan, '-', 'Color', c, 'LineWidth', 1.5);
% end
%
% yline(0, ':k', 'HandleVisibility', 'off');
% xlabel('Cranial caudal position (mm)');
% ylabel('Coherence difference');
% title('Significant coherence differences height sorted');
% legend(h, arrayfun(@(s) sprintf('Subj %d', s), 1:nSubj, 'UniformOutput', false), 'Location', 'bestoutside');
% set(gca, 'FontSize', 13)
% grid on;


    %% plot DICS results 2d
%     f=figure('Color','w'); hold on
%     diffCoh = source_perm.avgA.coh - source_perm.avgB.coh;
%     hLine = plot(sources_cent.pos(:,2), diffCoh, 'LineWidth', 2, 'Color', [0.2 0.6 0.8]);
%     hThresh = plot(sources_cent.pos(:,2), thr95, 'k--', 'LineWidth', 1.5);
%     aboveIdx = diffCoh > thr95;
%     xFill = [sources_cent.pos(aboveIdx,2); flipud(sources_cent.pos(aboveIdx,2))];
%     yFill = [thr95(aboveIdx); flipud(diffCoh(aboveIdx))];
%     hFill = fill(xFill, yFill, [0.9 0.3 0.3], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
%     xlabel('Cranial–caudal distance (mm)', 'FontWeight', 'bold', 'FontSize', 12)
%     ylabel('Contraction − Rest', 'FontWeight', 'bold', 'FontSize', 12)
%     grid on; box on
%     set(gca, 'FontSize', 12)
%     legend([hLine hThresh hFill], {'Difference', '95th percentile max', 'Significant'}, 'Location', 'best')
%     title(sprintf('%s coherence',sub))
