
function lay = make_layout_brainchans(ftdat,brainchan_labels)
%% make layout for brain
ftdat.grad.chanpos=ftdat.grad.coilpos;
ftdat.grad.chanori=ftdat.grad.coilori;

ftdat.grad.coordsys='ctf';

cfg=[];
cfg.grad=ftdat.grad;
cfg.channel=brainchan_labels(contains(brainchan_labels, 'Z'));
cfg.projection='orthographic';
cfg.viewpoint='auto';%topleft';%'anterior'; %topleft for brainspine
lay=ft_prepare_layout(cfg);

%ft_plot_layout(lay)



% % rotate deg
% theta = 270;  % Angle of rotation in degrees
% rotation_matrix = [cosd(theta), -sind(theta); sind(theta), cosd(theta)];
% 
% % Apply the rotation matrix to all electrode positions
% for i = 1:length(lay.pos)
%     % Rotate the electrode positions (x, y)
%     lay.pos(i, 1:2) = (rotation_matrix * lay.pos(i, 1:2)')';
% end
% 
% for i = 1:length(lay.mask{1})
%     lay.mask{1}(i, 1:2) = (rotation_matrix * lay.mask{1}(i, 1:2)')';
% end

ft_plot_layout(lay)
