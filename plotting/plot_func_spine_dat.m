function plot_func_spine_dat(mesh,src,func,grad,sensstl)

      inside=src.pos(src.inside,:);
      colmaptouse=brewermap([],'RdYlBu');
      ft_plot_mesh(mesh,'FaceColor', [0.8 0.8 1.0], 'facealpha', 0.2, 'EdgeColor','none',...
        'AmbientStrength', 0.15); hold on
      %camlight%('headlight');
  
      ft_plot_cloud(inside,func,'scalerad','no', 'cloudtype','surf','ncirc',1,'radius',5,'colormap',colormap(flipud(colmaptouse)))
      xlabel('X'); ylabel('Y')
      axis image
      

% Get the color limits for the first plot
cmin = min(func(:));
cmax = max(func(:));

% Determine tick positions
ticks = linspace(cmin, cmax, 4); % Change the number (here, 5) to position more ticks
tickLabels = arrayfun(@(x) sprintf('%.1f', x), ticks, 'UniformOutput', false);

% Apply colorbar settings
h = colorbar;
h.Ticks = ticks;
%h.TickLabels = tickLabels;



      hold on
      ft_plot_mesh(sensstl,'FaceColor', [0.7 0.7 1.0], 'facealpha', 1, 'EdgeColor','none',...
        'AmbientStrength', 0.15);
      view(7,1)
      camlight
      material('dull')
end



      %ft_plot_sens(grad,'coilshape', 'point','coilsize',18)

     % ft_plot_mesh(mesh, 'edgealpha',0.05, 'facealpha', 0.3); hold on




%       ft_plot_mesh(src.pos, 'vertexsize', t_value)%, 'vertexcolor', 'm');
% 
% 
%       ft_plot_cloud(pos, fun, 'mesh', anatomical,...
%         'radius', cfg.radius, 'rmin', cfg.rmin, 'scalerad', cfg.scalerad, ...
%         'ptsize', cfg.ptsize, 'ptdensity', cfg.ptdensity, 'ptgradient', cfg.ptgradient,...
%         'colorgrad', cfg.colorgrad, 'colormap', cfg.funcolormap, 'clim', [fcolmin fcolmax], ...
%         'unit', functional.unit, 'slice', cfg.slice, 'cloudtype', cfg.cloudtype, ...
%         'ori', cfg.ori, 'slicepos', cfg.slicepos, 'nslices', cfg.nslices, 'minspace', cfg.minspace,...
%         'intersectcolor', cfg.intersectcolor, 'intersectlinestyle', cfg.intersectlinestyle, ...
%         'intersectlinewidth', cfg.intersectlinewidth, 'ncirc', cfg.ncirc, ...
%         'scalealpha', cfg.scalealpha, 'facecolor', cfg.facecolor, 'edgecolor', cfg.edgecolor,...
%         'facealpha', cfg.facealpha, 'edgealpha', cfg.edgealpha, 'marker', cfg.marker,...
%         'vertexcolor', cfg.vertexcolor);
