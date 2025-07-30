
function lay = make_layout_neckchans(ftdat,sc_labels,ori)
%% make layout for sc
ftdat.grad.chanpos=ftdat.grad.coilpos;
ftdat.grad.chanori=ftdat.grad.coilori;

ftdat.grad.coordsys='ctf';


% need to remove bad chans from grad structrure
grad1=ftdat.grad; %copy
idx2keep=find(contains(grad1.label,sc_labels));

idx2keep = find(contains(grad1.label, sc_labels));

grad1.chanori = grad1.chanori(idx2keep, :);
grad1.chanpos = grad1.chanpos(idx2keep, :);
grad1.label = grad1.label(idx2keep);
grad1.chantype=grad1.chantype(idx2keep,:);
grad1.chanunit=grad1.chanunit(idx2keep,:);
grad1.coilori = grad1.coilori(idx2keep, :);
grad1.coilpos = grad1.coilpos(idx2keep, :);
grad1.tra=grad1.tra(idx2keep,idx2keep);

body.coordsys='ctf';



cfg=[];
cfg.grad=grad1;
cfg.headshape=body;
cfg.channel=grad1.label(contains(grad1.label, ori)); %radial
cfg.rotate=30
cfg.projection='orthographic';
%cfg.viewpoint= 'posterior';
lay=ft_prepare_layout(cfg);
ft_plot_layout(lay)


% % Define a transformation matrix (e.g., rotation around Z-axis)
% angle = deg2rad(-15); % adjust this angle until the projection looks correct
% R = [cos(angle), -sin(angle), 0, 0;
%      sin(angle),  cos(angle), 0, 0;
%      0,           0,          1, 0;
%      0,           0,          0, 1];
% 
% % Apply to grad structure
% grad1_rot = ft_transform_geometry(R, grad1);
% 
% % Optional: rotate body/headshape too
% body_rot = ft_transform_geometry(R, body);
% 
% % Then re-prepare layout
% cfg = [];
% cfg.grad = grad1_rot;
% cfg.headshape = body_rot;
% cfg.channel = grad1.label(contains(grad1.label, ori)); % radial
% cfg.projection = 'orthographic';
% cfg.viewpoint = 'posterior';
% lay = ft_prepare_layout(cfg);
% 
% ft_plot_layout(lay);


