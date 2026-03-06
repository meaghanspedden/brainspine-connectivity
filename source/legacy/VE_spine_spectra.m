% virtual electrode spinal cord (from prevalence analysis) 
clear all
close all
clc

addpath('C:\Users\mspedden\Documents\brainspineconnectivity\source')
addpath('C:\Users\mspedden\Documents\spm')
spm('defaults','EEG')

addpath('C:\Users\mspedden\Documents\fieldtrip')
ft_defaults

save_dir='C:\Users\mspedden\Documents\brainspine_save_newLF';
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
load('C:\Users\mspedden\Documents\brainspine_save_newLF\cluster_spineEMG_pos_newLF.mat')


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


   %% load and organise spinal cord leadfields
    %[Gx, Gy, Gz] = build_leadfield_matrices(fullfile(generic_dir,'cervical_realistic_brain_spine'), LFop);
    load('C:\Users\mspedden\Documents\bem_spine_fields\leadfield_cervical_realistic_bem_bem_.mat')

    nsourcepoints = size(leadfield_cord.pos,1);
    %nchannels     = size(Gx,2);
    spchanidx=find(grad_mm.coilpos(:,2) < 200); %indexed locally in grad struct (same indexing as LF)
    spchanlabs=grad_mm.label(spchanidx);

    %% clip leadfields to spinal cord channels only
%     Gx=Gx(:,spchanidx);
%     Gy=Gy(:,spchanidx);
%     Gz=Gz(:,spchanidx);
Lf=leadfield_cord;
Lf.label = leadfield_cord.label(spchanidx);
for i = 1:numel(leadfield_cord.leadfield)
    
    if ~isempty(leadfield_cord.leadfield{i})
        
        % leadfield{i} is [nChan × nOri]
        Lf.leadfield{i} = leadfield_cord.leadfield{i}(spchanidx, :);
        
    end
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

savename=sprintf('VE_spine_sub%s_forspectra_newLF', sub, 'VE');
save(fullfile(save_dir,savename), 'VE')


end