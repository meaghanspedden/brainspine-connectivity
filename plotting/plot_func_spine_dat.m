function plot_func_spine_dat(mesh,src,func,grad,senscol)

      inside=src.pos(src.inside,:);
      colmaptouse=brewermap([],'YlGnBu');
      ft_plot_mesh(mesh, 'edgealpha',0.05, 'facealpha', 0.3); hold on
      ft_plot_cloud(inside,func,'scalerad','no', 'cloudtype','surf','ncirc',1,'radius',5,'colormap',colormap(flipud(colmaptouse)))
      xlabel('X'); ylabel('Y')
      axis image

      hold on
      ft_plot_sens(grad,'coilshape','sphere','coilsize',10)
      view(7,1)
      colorbar
end








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
