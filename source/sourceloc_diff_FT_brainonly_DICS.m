
%% script that runs DICS in FT to image brain coherence with ref chan (EMG) 

clear all
close all
clc

addpath('C:\Users\mspedden\Documents\brainspineconnectivity\source')
addpath('C:\Users\mspedden\Documents\spm')
spm('defaults','EEG')

addpath('C:\Users\mspedden\Documents\fieldtrip')
ft_defaults

save_dir='C:\Users\mspedden\Documents\brainspine_save';
if ~exist(save_dir,'dir')
    mkdir(save_dir)
end

% exclude two subejcts because one had very small head (very far away from
% sensors) and other couldnt close headcast

subs = {'OP00212','OP00213','OP00215', 'OP00219', ...
 'OP00225', 'OP00221', 'OP00224'};

generic_dir = 'C:\Users\mspedden\Documents\new_leadfields_and_geom'; %where I have saved folder with brainspione leadfields and geoms (meshes)
geomfile = fullfile(generic_dir, 'geometries_cervical_realistic.mat');
rng(1)

HFC=1;
rectify=1;
mult_comp_corr=1;
fband=[10 35];

subjResults=struct();

for ss=1:length(subs)

    sub=subs{ss};

    if strcmp(sub, 'OP00224')
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

    if rectify
        cfg=[];
        cfg.rectify='yes';
        cfg.channel='EXG1';
        ftdatr=ft_preprocessing(cfg,ftdat);

        for k=1:length(ftdat.trial)
            ftdat.trial{k}(end,:)=ftdatr.trial{k}; %ftdat has rectified emg
        end

    end

    brainlabs=[];
    for f = 1:length(grad_mm.label)
        if grad_mm.chanpos(f,2) > 200
            brainlabs = [brainlabs grad_mm.label(f)];
        end
    end

    %% new brain only dataset and grad

brainidx = find(grad_mm.chanpos(:,2) > 200);
braingrad = grad_mm;
braingrad.chanori = grad_mm.chanori(brainidx, :);
braingrad.chanpos = grad_mm.chanpos(brainidx, :);
braingrad.chantype = grad_mm.chantype(1:length(brainidx));
braingrad.chanunit = grad_mm.chanunit(1:length(brainidx));
braingrad.coilori = grad_mm.coilori(brainidx, :);
braingrad.coilpos = grad_mm.coilpos(brainidx, :);
braingrad.label = grad_mm.label(brainidx);
braingrad.tra = grad_mm.tra(brainidx, brainidx);

%figure; ft_plot_sens(braingrad)


%use all channels but only brain grad

braindat=ftdat;
mesh_brain.unit='mm';

%single shell volume conductor
cfg = [];
cfg.method = 'singleshell';
vol = ft_prepare_headmodel(cfg, mesh_brain);

% calculate forward
cfg                     = [];
cfg.sourcemodel         = sources_brain;
cfg.headmodel           = vol;
cfg.grad                = braingrad;
cfg.reducerank          = 'no';

LF = ft_prepare_leadfield(cfg);


    %% beamforming----------------------------
    cfg=[];
    cfg.output     = 'powandcsd';
    cfg.method     = 'mtmfft';
    cfg.foilim     = fband;
    cfg.tapsmofrq  = 1;
    cfg.keeptrials='yes';
    freqdat_tr=ft_freqanalysis(cfg,braindat); %trial wise freq

    cfg=[];
    cfg.avgoverfreq='yes';
    freqdat_tr=ft_selectdata(cfg,freqdat_tr);%need to average over frequency as input to permutation

    %% separate conditions

    statidx=find(ftdat.trialinfo==1);
    restidx=find(ftdat.trialinfo==2);

    [nTrials,~] = min([length(statidx) length(restidx)]);

    cfg=[];
    cfg.trials=statidx(1:nTrials);
    statdat=ft_selectdata(cfg,freqdat_tr); %trial wise data pr condition

    cfg.trials=restidx(1:nTrials);
    restdat=ft_selectdata(cfg,freqdat_tr);

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

%filter from both conditions
    cfg=[];
    cfg.grid = sources_brain;
    cfg.headmodel=vol;
    cfg.sourcemodel.leadfield=LF;
    cfg.dics.keepfilter='yes';
    cfg.dics.lambda=10;
    cfg.method = 'dics';
    cfg.refchan='EXG1';

    coh_source = ft_sourceanalysis(cfg,freqdat);

    cfg=[];
    cfg.grid = sources_brain;
    cfg.headmodel=vol;
    cfg.sourcemodel.leadfield=LF;
    cfg.dics.filter=coh_source.avg.filter;
    cfg.dics.lambda=10;
    cfg.method = 'dics';
    cfg.refchan='EXG1';

    %source_stat = ft_sourceanalysis(cfg,statdat); %struct per condition with single trials
    % source_rest = ft_sourceanalysis(cfg,restdat);
    
    cfg.permutation = 'yes';
    cfg.numpermutation=500;
    source_perm = ft_sourceanalysis(cfg, statdat, restdat);
   
    nsourcepoints=length(source_perm.inside);
    nPerm = numel(source_perm.trialA);
    cohDiff_perm = zeros(nsourcepoints, nPerm);

    for i = 1:nPerm
        cohA_perm = source_perm.trialA(i).coh;
        cohB_perm = source_perm.trialB(i).coh;
        cohDiff_perm(:, i) = cohA_perm - cohB_perm;
    end

    maxPerm = max(cohDiff_perm, [], 1); % max over sources
    
    if mult_comp_corr
        thr95 = prctile(maxPerm, 95,2);   % family-wise threshold (dim 1x1)
    else
        thr95=(ones(1,nsourcepoints)*thr95)'; %otherwise get one per source point
    end
   
    coh_diff=source_perm.avgA.coh-source_perm.avgB.coh;
    source_diff=coh_source;
    source_diff.avg.coh=coh_diff; %raw coh difference

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
    brain_int = ft_sourceinterpolate(cfg, source_p, mesh_brain);

    
%     permMean = mean(cohDiff_perm, 2);
%     permStd  = std(cohDiff_perm, 0, 2);
%     
%     z_coh = (coh_diff - permMean) ./ permStd;
%     permStd(permStd == 0) = NaN;
% 
%     source_z = coh_source;        %source structure
%     source_z.avg.coh = z_coh; 

%     cfg = []; %interpolate raw coherence diff
%     cfg.parameter = {'coh'};
%     brain_int=ft_sourceinterpolate(cfg,source_diff, mesh_brain);
% 
%     cfg = []; %interpolate z score
%     cfg.parameter = {'coh'};
%     brain_int_z=ft_sourceinterpolate(cfg,source_z, mesh_brain);
    
    %% add a mask for stat significance
    source_mask=coh_source; %copy
    source_mask.avg.pow = coh_diff > thr95;

%     %interpolate the mask
%     cfg = [];
%     cfg.parameter = 'pow';
%     cfg.interpmethod = 'nearest';
%     source_mask_int = ft_sourceinterpolate(cfg, source_mask, mesh_brain);
%     brain_int.mask=source_mask_int.pow;
%     
%     brain_int_z.mask=source_mask_int.pow;

source_mask = coh_source;
source_mask.avg.coh = double(mask);  % logical → double for FT
cfg = [];
cfg.parameter = 'coh';
cfg.interpmethod = 'nearest';
source_mask_int = ft_sourceinterpolate(cfg, source_mask, mesh_brain);
%     
% brain_int.mask = source_mask_int.coh;  % same as before


% mask = source_mask.avg.pow;  % logical array same size as coh_diff
% 
% % Find the max only where mask == true (not using this pt)
% coh_diff_masked = coh_diff;
% coh_diff_masked(~mask) = -Inf;  
% [max_val, max_idx] = max(coh_diff_masked(:));
% 
% % Now  map back to coordinates
% max_pos = sources_brain(max_idx, :);

ncol = 256;
addpath('C:\Users\mspedden\Documents\fieldtrip\external\matplotlib\') 
brain_color = [0.92 0.92 0.92];   % first color = background
hotmap = flipud(magma(ncol-1));       

cmap = [brain_color; hotmap];
%     
    f= figure;
    cfg = [];
    cfg.figure='gcf';
    cfg.method = 'surface';
    cfg.funparameter = 'coh';
    cfg.surfalpha=.3;
    cfg.funcolormap = cmap;
    cfg.funcolorlim  = [0 max(brain_int.coh(:))];  % 0 maps to brain_color
    cfg.projmethod = 'nearest';
    cfg.surffile = 'mesh_brain';
    ft_sourceplot(cfg, brain_int);
    view(176,-10)
    camlight
    ax = gca;                  
    ax.FontSize = 14;
hpatch = findobj(gcf, 'Type', 'patch');
set(hpatch, 'FaceAlpha',0.9)

    subjResults(ss).coh_diff = coh_diff; %orig brain space (not interpolated)
    subjResults(ss).sig_mask = source_mask_int.coh; % binary 0/1 map
    subjResults(ss).pos = source_mask_int.pos; 
    subjResults(ss).inside = source_mask_int.inside;
    subjResults(ss).thr95=thr95;


end


%cd C:\Users\mspedden\Documents\brainspine_save
load('groupRes_brain_DICS.mat')

save(fullfile(save_dir,'groupRes_brain_DICS'), 'subjResults')

nSubs = length(subjResults);

% cat all binary significance masks
all_masks = cat(2, subjResults(:).sig_mask);
sig_pos = false(nSubs,1);

%number of subjects with at least one significant source anywhere
for ss = 1:nSubs
    diffCoh = subjResults(ss).coh_diff;   % source_perm.avgA.coh - avgB.coh
    thr95    = subjResults(ss).thr95; % 95th percentile from permutation

    if any(diffCoh > thr95)
        sig_pos(ss) = true;
    end
end

fprintf('%g out of %g subjects show sig coherence in brain\n', sum(sig_pos), nSubs)

% binomial p and CI
x = sum(sig_pos);
n = nSubs;
alpha=0.05;
phat = x/n;
lower = betainv(alpha/2, x,     n-x+1);
upper = betainv(1-alpha/2, x+1, n-x);

ci = [lower upper];
p = 1 - binocdf(x-1, n, alpha);


%% prevalence sig res
group_prevalence = mean(all_masks, 2);

%  group source structure for plotting
group_source = subjResults(1); % use one subject as template
group_source.pow = group_prevalence;
group_source = rmfield(group_source, {'coh_diff', 'sig_mask'}); % clean up

group_ft = []; %combine this with code above
group_ft.pos = group_source.pos;
group_ft.inside = group_source.inside;
group_ft.pow = group_prevalence;


%% Interpolate group map onto the brain mesh
cfg = [];
cfg.parameter = 'pow';
cfg.interpmethod = 'nearest';
group_int = ft_sourceinterpolate(cfg, group_ft, mesh_brain);

threshold = 0.3;
group_prevalence_masked = group_prevalence;
group_prevalence_masked(group_prevalence < threshold) =0;  % or 0

group_ft.pow    = group_prevalence_masked;

cfg = [];
cfg.parameter = 'pow';
cfg.interpmethod = 'nearest';
group_int = ft_sourceinterpolate(cfg, group_ft, mesh_brain);

% group_int.pow_thresh = group_int.pow;          % copy original
% group_int.pow_thresh(group_int.pow < threshold) = 0;   % subthreshold = NaN

%% Plot group prevalence map
figure;
cfg = [];
cfg.method = 'surface';
cfg.funparameter = 'pow';
cfg.funcolorlim =  [threshold max(group_int.pow)];
cfg.funcolormap = cmap;
cfg.projmethod = 'nearest';
cfg.surffile = mesh_brain;
ft_sourceplot(cfg, group_int);
view(176, -10);
camlight;
ax = gca;                  % get current axes
ax.FontSize = 14;
hpatch = findobj(gcf, 'Type', 'patch');
set(hpatch, 'FaceAlpha',0.9)


%% binomial p and CI
x = sum(sig_pos);
n = nSubjects;
alpha=0.05;
phat = x/n;
lower = betainv(alpha/2, x,     n-x+1);
upper = betainv(1-alpha/2, x+1, n-x);

ci = [lower upper];
p = 1 - binocdf(x-1, n, alpha);

