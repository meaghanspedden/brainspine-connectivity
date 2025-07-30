
%evaluate lead fields

sub='OP00215';
subject_dir = fullfile('D:\MSST001', [sub '_merged']); 

which_ori=3;


rootlf   = fullfile(subject_dir, 'leadfields_nobrainchans.mat');
load(fullfile(subject_dir,'grad1.mat'));
sensor_data=grad1;

lf=load(rootlf);

leadfields=lf.leadfield.leadfield;
nsources=length(leadfields);

lfmat=[];

for k=1:nsources
    lfmat(k,:)=leadfields{k}(:,which_ori);
end

p = eye(size(lfmat,2)) - (sensor_data.chanori*pinv(sensor_data.chanori));
lfmat=lfmat';
final_out = p*lfmat;


supp=10*log10(var(final_out)./var(lfmat)); %dB each source is suppressed by
figure
histogram(supp)
xlabel('dB suppression')
ylabel('count (sourcepoints)')
title(sprintf('hfc effects brainspine leadfields ori %g',which_ori))

%(should be neg or zero)

%compare to shielding factor 

