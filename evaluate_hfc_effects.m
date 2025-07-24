
%evaluate lead fields

sub='OP00215';
subject_dir = fullfile('D:\MSST001', [sub '_static']); 



filename_static = fullfile('D:\MSST001', ...
    ['sub-' sub], ...
    'ses-001', ...
    'meg', ...
    ['pstaticoe1000mspddfflo45hi45hfcstatic_001_array1.mat']);

rootlf   = fullfile(subject_dir, 'leadfields.mat');
geomfile = fullfile(subject_dir, 'geoms.mat');
castforward = load(geomfile);
lf=load(rootlf);

D=spm_eeg_load(filename_static);
sensor_data=D.sensors('MEG');

leadfields=lf.leadfield.leadfield;
nsources=length(leadfields);
lfmat=[];
which_ori=1;
for k=1:nsources

lfmat(k,:)=leadfields{k}(:,which_ori);
end

p = eye(nsources) - (sensor_data.chanori*pinv(sensor_data.chanori));
final_out = p*lf_mat;




corr_mat = [];
for i=1:size(final_out,2)
    val = corrcoef(final_out(:,i), final_lf(:,i));
    corr_mat = [corr_mat val(2,1)];
end




    % Loop over sources
    figure;
    t = tiledlayout(1,3, 'Padding', 'compact', 'TileSpacing', 'compact');
    for k=1:3
        ax = nexttile(t);
        src_val = dip_inds_x(src_interest(k))+j-1;
        plot_topoplot_xz(ax, back_pos, final_out(endsWith(grad_back.label, '-TAN'), src_val), cmap)
        cl = clim;
        clim([-1 1]*max(abs(cl)));
        hold on;
        scatter(ax, transformedmesh(:,1)/1000, transformedmesh(:,3)/1000)
        title(sprintf('%s Source', source{k}));
    
    end
    title(t, sprintf('%s Orientated Sources - P*Lead fields', orientations{j}))

