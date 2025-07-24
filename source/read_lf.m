function [L,gainmatnames,lflabels]=read_lf(rootlf,Ncordpoints,castforward)
 %% lead fields for each source orientation from Maike format
for f=1:Ncordpoints,
    a1=load([rootlf sprintf('bem_%d_ori1.mat',f)]); %Ls_x 1st 64 channels, Ls_y 2nd 64 channels
    a2=load([rootlf sprintf('bem_%d_ori2.mat',f)])
    a3=load([rootlf sprintf('bem_%d_ori3.mat',f)])
    L(1,f,:)=[a1.Ls_x;a1.Ls_y;a1.Ls_z];
    L(2,f,:)=[a2.Ls_x;a2.Ls_y;a2.Ls_z];
    L(3,f,:)=[a3.Ls_x;a3.Ls_y;a3.Ls_z];
end;
chanpos_cast=castforward.coils_3axis.chanpos;
lflabels=castforward.coils_3axis.label;
label=lflabels;
for f=1:3,
    gainmatnames{f} = ['SPMgainmatrix_' sprintf('cord%d_grb',f)  '.mat'];
    G=squeeze(L(f,:,:));
    save([rootlf filesep gainmatnames{f}], 'G', 'label', spm_get_defaults('mat.format'));
end;


L=[]; %% lead fields for each source orientation
for f=1:Ncordpoints,
    a1=load([rootlf sprintf('bem_%d_ori1.mat',f)]); %Ls_x 1st 64 channels, Ls_y 2nd 64 channels
    a2=load([rootlf sprintf('bem_%d_ori2.mat',f)])
    a3=load([rootlf sprintf('bem_%d_ori3.mat',f)])
    L(f,:,1)=[a1.Ls_x;a1.Ls_y;a1.Ls_z];
    L(f,:,2)=[a2.Ls_x;a2.Ls_y;a2.Ls_z];
    L(f,:,3)=[a3.Ls_x;a3.Ls_y;a3.Ls_z];
end;
chanpos_cast=castforward.coils_3axis.chanpos;
lflabels=castforward.coils_3axis.label;
label=lflabels;
for f=1:3,
    gainmatname = ['SPMgainmatrix_' sprintf('cord%d_grb',f)  '.mat'];
    G=squeeze(L(:,:,f));
    save([rootlf filesep gainmatname], 'G', 'label', spm_get_defaults('mat.format'));
end;