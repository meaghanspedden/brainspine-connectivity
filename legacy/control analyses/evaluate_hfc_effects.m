
%evaluate lead fields

sub='OP00212';
analysis='merged';
subject_dir = fullfile('D:\MSST001', [sub '_merged']); 

which_ori=3;

[Gx, Gy, Gz] = build_leadfield_matrices('D:\new_leadfields_and_geom\cervical_cord_spineonly');


filename = fullfile('D:\MSST001', ...
    ['sub-' sub], ...
    'ses-001', ...
    'meg', ...
    ['p' analysis 'oe1000mspddfflo45hi45hfcstatic_001_array1.mat']);

D=spm_eeg_load(filename);

nsources=size(Gx,1);

sensor_data=D.sensors('MEG');

sens_pos=sensor_data.coilpos(:,2);

spineind=[];
for k=1:length(sens_pos)
        if sens_pos(k) < 200
            spineind = [spineind k];
        end
end



if which_ori==1
    lfmat=Gx;

elseif which_ori==2
    lfmat=Gy;
else
    lfmat=Gz;
end


lfmat=lfmat(:,spineind);
chanori=sensor_data.coilori(spineind,:);

load("M.mat")
p=M(spineind,spineind);

%lfmat needs to be nsources x nchannels
%p = eye(size(lfmat,2)) - (chanori*pinv(chanori));
lfmat=lfmat';
final_out = p*lfmat;


supp=10*log10(var(final_out)./var(lfmat)); %dB each source is suppressed by
figure
histogram(supp)
xlabel('dB suppression')
ylabel('count (sourcepoints)')
title(sprintf('hfc effects spinal cord leadfields ori %g',which_ori))

%(should be neg or zero)

%compare to shielding factor 

