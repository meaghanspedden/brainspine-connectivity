function plotOptOri(Jv,src,subject)

figure;
ft_plot_mesh(subject,'EdgeAlpha',0.2,'FaceAlpha',0.2); hold on
%ft_plot_mesh(cyl,'Ed')
source_id = find(src.inside);
orix=ones(length(source_id),1)*Jv(1);
oriy=ones(length(source_id),1)*Jv(2);
oriz=ones(length(source_id),1)*Jv(3);
quiver3(src.pos(source_id,1),src.pos(source_id,2),src.pos(source_id,3),orix, oriy,oriz,'LineWidth',2);
%scatter3(src.pos(source_id,1),src.pos(source_id,2),src.pos(source_id,3),'r*');
axis image
xlabel('X'); ylabel('Y'); zlabel('Z')
view(0,0)
title(' Opt orientation')
