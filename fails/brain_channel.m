%% brain channel analysis

clear all
close all

sub='OP00215';

addpath('D:\msg_source_recon\')

save_dir = fullfile('D:\MSST001', [sub '_contrast']);
if ~exist(save_dir, 'dir')
    mkdir(save_dir)
end
subject_dir = fullfile('D:\MSST001', [sub '_static']); %rest and static will be same here

filename_rest = fullfile('D:\MSST001', ...
    ['sub-' sub], ...
    'ses-001', ...
    'meg', ...
    ['prestoe1000mspdfflo45hi45hfcrest_001_array1.mat']);

filename_static = fullfile('D:\MSST001', ...
    ['sub-' sub], ...
    'ses-001', ...
    'meg', ...
    ['pstaticoe1000mspdfflo45hi45hfcstatic_001_array1.mat']);

geomfile = fullfile('D:\MSST001', [sub '_static'], 'geoms.mat');


%% fixed across participants
brainlabs={'41', '48', '34', '43', '36', '10', '12',...
    '44', '45', '33', '15', '39', '47'};

%brainlabs=labelToIndex(brainchan_labels);


% LOAD  DATA--------------------------------------------------------
Drest = spm_eeg_load(filename_rest);
Dstatic = spm_eeg_load(filename_static);
load(geomfile)

grad=Drest.sensors('MEG');

[powTrialS, chanlabs2use] =getPow(Dstatic,brainlabs);
[powTrialR, ~] =getPow(Drest,brainlabs);

gradbrain=find(contains(grad.label,chanlabs2use));


gradb=struct();
gradb.chanori=grad.chanori(gradbrain,:);
gradb.chanpos=grad.chanpos(gradbrain,:);
gradb.chantype=grad.chantype(gradbrain,:);
gradb.chanunit=grad.chanunit(gradbrain,:);
gradb.coilori=grad.coilori(gradbrain,:);
gradb.coilpos=grad.coilpos(gradbrain,:);
gradb.label=grad.label(gradbrain,:);
gradb.tra=grad.tra(gradbrain,gradbrain);
gradb.unit=grad.unit;
gradb.balance=grad.balance;

%% t test rest static for each z channel
ntrials=min([length(powTrialS) length(powTrialR)]);
tvals=[]; pvals=[];

for c=1:size(powTrialS,2) %loop through channels
    statdat=log(powTrialS(:,c));
    restdat=log(powTrialR(:,c));
    [h,p,ci,stats] = ttest2(statdat,restdat);
    tvals(c)=stats.tstat^2;
    pvals(c)=p;

end


[pkpos, peakind_pos] = max(tvals);  
pkData=powTrialS(:,peakind_pos);

brainlabmax=chanlabs2use(peakind_pos);
[~, ~, ~, adj_p]=fdr_bh(pvals,0.05,'dep','no');
pSig=adj_p<0.05;


ft_plot_mesh(mesh_torso,'FaceColor', [0.8 0.8 1.0], 'facealpha', 0.2, 'EdgeColor','none',...
    'AmbientStrength', 0.15); hold on
ft_plot_cloud(gradb.chanpos,tvals,'scalerad','no', 'cloudtype','surf','radius',3.5)
xlabel('X'); ylabel('Y')
axis image

plot3(gradb.chanpos(peakind_pos,1), gradb.chanpos(peakind_pos,2), gradb.chanpos(peakind_pos,3), 'ro', 'MarkerSize',12)
view(-0.8428,-1.8897)
colorbar

 % Save results
 bffilename = fullfile(save_dir, ['chandata_brain_', sub]);

 save(bffilename, 'pkData', 'powTrialS', 'tvals','pSig','chanlabs2use', 'brainlabmax');















% %make into 3d volume for nifti printing
% nX = length(unique(allhp(:,1)));
% nY = length(unique(allhp(:,2)));
% nZ = length(unique(allhp(:,3)));
%
% %map grid positions to volume indices
% xvals = unique(allhp(:,1));
% yvals = unique(allhp(:,2));
% zvals = unique(allhp(:,3));
%
% voxsize = [xvals(2)-xvals(1), yvals(2)-yvals(1), zvals(2)-zvals(1)];
% origin  = [xvals(1), yvals(1), zvals(1)];
%
% affine = [voxsize(1), 0,          0,          origin(1);
%          0,          voxsize(2), 0,          origin(2);
%          0,          0,          voxsize(3), origin(3);
%          0,          0,          0,          1];
%
% V = struct();
% V.dim = [nX nY nZ];
% V.dt = [16 0];  % float32
% V.mat = affine;
% V.descrip = 'Functional source data';
%
% for im=1:5 %length(source_brain_tr.trial)
% % Step 3: For each inside point, get indices in volume
% vol3D = nan(nX, nY, nZ);
%
% for i = 1:size(allhp,1)
%     if isinside(i)
%         xidx = find(xvals == allhp(i,1));
%         yidx = find(yvals == allhp(i,2));
%         zidx = find(zvals == allhp(i,3));
%         vol3D(xidx, yidx, zidx) = source_brain_tr.trial(im).pow(i);
%     end
% end
%
%
% V.fname = sprintf('trial_%g_%s.nii',im,analysis);
% cd(save_dir)
% spm_write_vol(V, vol3D);
%
%
% end



% tic
% brain_source=[]; %loop through src pts and trials
% for s=1:length(brainsrc.pos)
%
%     for t=1:length(ftdat.trial)
%         bf_filt = source.avg.filter{s};
%         brain_source(s,:,t)= bf_filt * ftdat.trial{t}(brainchanidx,:);
%     end
%
% end
% toc





   


