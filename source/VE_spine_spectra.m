% virtual electrode spinal cord (from prevalence analysis) 
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

%subs={'OP00212'};
% exclude two subejcts because one had very small head (very far away from
% sensors) and other couldnt close headcast

subs = {'OP00212','OP00213','OP00215', 'OP00219', ...
 'OP00225', 'OP00221', 'OP00224'}; %have brain 

generic_dir = 'C:\Users\mspedden\Documents\new_leadfields_and_geom'; %where I have saved folder with brainspione leadfields and geoms (meshes)
geomfile = fullfile(generic_dir, 'geometries_cervical_realistic.mat');

LFop='spine';
HFC=1;
rectify=1;
fband=[10 35];

%load structure with position
load('C:\Users\mspedden\Documents\brainspine_save\cluster_spineEMG_pos.mat')


% cmap = magma(256);          
% pointColor = cmap(90,:);
% meshColor = cmap(40,:);   % dark purple/black end of magma
% outlineColor=cmap(240,:);
% facecol=[0.72 0.70 0.61];

% figure
% ft_plot_mesh(mesh_wm, ...
%     'facecolor', facecol, ...   % FT default: light gray
%     'facealpha', 1, ...               % opaque by default
%     'edgecolor', 'none');             % no edges
% hold on
% view(-250, -1)
% lighting gouraud
% camlight


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

    %remove bad channels
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


    [Gx, Gy, Gz] = build_leadfield_matrices(fullfile(generic_dir,'cervical_realistic_brain_spine'), LFop);

    nsourcepoints = size(Gx,1);
    nchannels     = size(Gx,2);
    spchanidx=find(grad_mm.coilpos(:,2) < 200); %indexed in grad struct
    spchanlabs=grad_mm.label(spchanidx);


 %% clip leadfields to spinal cord channels only
    Gx=Gx(:,spchanidx);
    Gy=Gy(:,spchanidx);
    Gz=Gz(:,spchanidx);

    Lf.pos    = sources_cent.pos;     % nsourcepoints x 3
    Lf.inside = sources_cent.inside;     % all points inside
    Lf.unit   = 'mm';
    Lf.label  = grad_mm.label(spchanidx);   % nchannels x 1 cell
    Lf.leadfielddimord = '{pos}_chan_ori';

    Lf.leadfield = cell(1,nsourcepoints);

    for k = 1:nsourcepoints
        % Combine X/Y/Z components
        Lf.leadfield{k} = [Gx(k,:)' Gy(k,:)' Gz(k,:)']; % nchannels x 3
    end

    % 2. dummy head model for input config only
    cfg                     = [];
    cfg.method              = 'infinite';
    cfg.siunits=1;
    cfg.grad=grad_mm;
    cfg.conductivity = 1;

    dummyvol = ft_prepare_headmodel(cfg,mesh_torso);
   

%% get pos and index for point with max coherence with EMG
if ss==1
f = figure('Color', 'w'); % white background
spineColor = [0.7 0.7 0.7];

ft_plot_mesh(mesh_wm, ...
    'facecolor', spineColor, ...
    'facealpha', 0.3,...
    'edgecolor', 'none'); 
hold on
highlightColor = [1 0.4 0.1];
plot3(ROIpos(:,1), ROIpos(:,2), ROIpos(:,3), 'o', ...
    'MarkerSize', 10, ...
    'MarkerEdgeColor', [0.9 0.3 0], ...
    'MarkerFaceColor', highlightColor, ...
    'LineWidth', 2);
[sx, sy, sz] = sphere(40);
r = .5;
view(90,18)
material dull

end

%% lcmv

cfg                   = [];
cfg.covariance        = 'yes';
tlock = ft_timelockanalysis(cfg, ftdat);

cfg                     = [];
cfg.method              = 'lcmv';
cfg.headmodel           = dummyvol;
cfg.sourcemodel.leadfield     = Lf;
cfg.grid     = sources_cent;
cfg.unit                = Lf.unit;
cfg.lcmv.keepfilter     = 'yes';
source_idx = ft_sourceanalysis(cfg, tlock);

roi_center = mean(ROIpos, 1);
d = sqrt(sum((ROIpos - roi_center).^2, 2));
R = max(d);

d_all = sqrt(sum((source_idx.pos - roi_center).^2, 2));
roi_idx_radius = find(d_all <= R & source_idx.inside);

cfg = [];
cfg.pos          = roi_center;   
cfg.radius       = R;           
cfg.method       = 'svd';
cfg.numcomponent = 1;
VE = ft_virtualchannel(cfg, ftdat, source_idx);

savename=sprintf('VE_spine_sub%s_forspectra', sub, 'VE');
save(fullfile(save_dir,savename), 'VE')
%%


% beamformer=source_idx.avg.filter{1};
% 
% chansel  = tlock.label(1:end-1); % MEG sensors
% chanindx = 1:length(chansel);    %hard coded    
% 
% coh_data = [];
% coh_data.label = {'coh_x', 'coh_y', 'coh_z'};
% coh_data.time = braindat.time;
% 
% %this gives for each orientation
% for i=1:length(braindat.trial)
%   coh_data.trial{i} = beamformer * braindat.trial{i}(chanindx,:);
% end
% 
% timeseries = cat(2, coh_data.trial{:});
% [u1, s1, v1] = svd(timeseries, 'econ');
% 
% virtualchannelbrain = [];
% virtualchannelbrain.label = {'motor'};
% virtualchannelbrain.time = braindat.time;
% 
% for k = 1:length(braindat.trial)
%   virtualchannelbrain.trial{k}(1,:) = u1(:,1)' * beamformer * braindat.trial{k}(chanindx,:);
% end
% 
% savename=sprintf('brain_VC_%s',sub);
% save(fullfile(save_dir,savename), 'virtualchannelbrain')
% error('stop')
% 
% % cfg = [];
% % cfg.channel = 'EXG1';
% % emgdata = ft_selectdata(cfg, braindat);
% % 
% % cfg = [];
% % combineddata = ft_appenddata(cfg, virtualchannelbrain,virtualchanneldata);
% % 
% % cfg            = [];
% % cfg.output     = 'fourier';
% % cfg.method     = 'mtmfft';
% % cfg.foilim     = [5 45];
% % cfg.tapsmofrq  = 3;
% % cfg.keeptrials = 'yes';
% % freq = ft_freqanalysis(cfg, combineddata);
% % 
% % cfg = [];
% % cfg.method = 'coh';
% % coherence = ft_connectivityanalysis(cfg, freq);
% % figure; plot(coherence.freq,(squeeze(coherence.cohspctrm(1,:,:))'))
% % ylim([0 0.15])
% 
% % savename=sprintf('sub%s_VE_brain_EMG',sub);
% % save(fullfile(save_dir,savename),"combineddata")

end