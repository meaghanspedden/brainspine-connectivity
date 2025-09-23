function [Lspine, Lbrain, gainmatnames, lflabels] = read_lf(rootlf, Nspine, Nbrain, castforward)
%% READ_LF - Load leadfields for spine and brain separately

nsensors = numel(castforward.coils_3axis.label)/3;
nchannels = nsensors * 3; 

Lspine = zeros(3, Nspine, nchannels);
Lbrain = zeros(3, Nbrain, nchannels);


%% --------- Spine Leadfields ---------
for idx = 1:Nspine
    for ori = 1:3
        S = load(fullfile(rootlf, sprintf('spine_%d_ori%d.mat', idx, ori)));
        Lspine(ori, idx, :) = [S.Ls_x; S.Ls_y; S.Ls_z];
    end
end


%% --------- Brain Leadfields ---------
for idx = 1:Nbrain
    for ori = 1:3
        S = load(fullfile(rootlf, sprintf('brain_%d_ori%d.mat', idx, ori)));
        Lbrain(ori, idx, :) = [S.Ls_x; S.Ls_y; S.Ls_z];
    end
end


%% --------- Channel Info ---------
lflabels = castforward.coils_3axis.label;

%% --------- Save Gain Matrices ---------
gainmatnames = {};
for ori = 1:3
    % Spine
    Gspine = squeeze(Lspine(ori,:,:));
    gname_spine = sprintf('SPMgainmatrix_spine_ori%d.mat', ori);
    save(fullfile(rootlf, gname_spine), 'Gspine', 'lflabels', ...
        spm_get_defaults('mat.format'));
    gainmatnames{end+1} = gname_spine;

    % Brain
    Gbrain = squeeze(Lbrain(ori,:,:));
    gname_brain = sprintf('SPMgainmatrix_brain_ori%d.mat', ori);
    save(fullfile(rootlf, gname_brain), 'Gbrain', 'lflabels', ...
        spm_get_defaults('mat.format'));
    gainmatnames{end+1} = gname_brain;
end
