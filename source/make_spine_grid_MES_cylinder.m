function srcOut = make_spine_grid_MES_cylinder(subject,cyl,res)


% make grid in cylinder
minpos = min(cyl.pos);
maxpos = max(cyl.pos);


cfg = [];
cfg.method = 'basedongrid';
cfg.xgrid = minpos(1):res:maxpos(1);
cfg.ygrid = minpos(2):res:maxpos(2);
cfg.zgrid = minpos(3):res:maxpos(3);
cfg.unit = 'mm';
warning('source grid unit hard coded to mm')
src = ft_prepare_sourcemodel(cfg);

for ii = 1:length(src.pos) %Georges code here
    
    tmp = src.pos(ii,:);
    sa = solid_angle(cyl.pos-tmp,cyl.tri);
    
    if abs(abs(sum(sa)) - 4*pi) < 1e-6
        src.inside(ii) = 1;
    else
        src.inside(ii) = 0;
    end
end

%easier just to deal with inside points
src.pos=src.pos(src.inside,:); 
src.inside=true(length(src.pos),1);

%% Plot
% figure;
% ft_plot_mesh(subject,'EdgeAlpha',0.2,'FaceAlpha',0.2); hold on
% source_id = find(src.inside);
% 
% scatter3(src.pos(source_id,1),src.pos(source_id,2),src.pos(source_id,3),'r*');
% axis image
% xlabel('X'); ylabel('Y'); zlabel('Z')


fprintf('Source grid has %g points, inside only\n',length(src.pos))
srcOut=src;



