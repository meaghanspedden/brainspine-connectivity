figure;
    ft_plot_mesh(mesh_bone, 'facecolor', [0.9 0.85 0.7], 'facealpha', 0.3, 'edgecolor', 'none');
hold on
plot3(sources_cent.pos(:,1), sources_cent.pos(:,2), sources_cent.pos(:,3),'ko')

nsourcepoints = size(sources_cent.pos,1);

for s = 1:nsourcepoints
    text(sources_cent.pos(s,1), ...
         sources_cent.pos(s,2), ...
         sources_cent.pos(s,3) + 2, ...  % small z-offset so text floats above
         num2str(s), ...
         'FontSize', 10, ...
         'Color', 'k', ...
         'HorizontalAlignment','center');
end

    view(-250, -1);


figure; plot(xpos,invp,'-o'); hold on
yline(invpthr, 'k--')
for s = 1:length(xpos)
    text(xpos(s), ...
         invp(s) + 0.02*range(invp), ...  % small vertical offset
         num2str(s), ...
         'HorizontalAlignment','center', ...
         'FontSize', 10);
end