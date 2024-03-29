Lfnorm=sqrt(sum(abs(Lf).^2, 1));
%our orientation only
Lfnorm_ori=Lfnorm(xind).*Jv(1)+Lfnorm(yind).*Jv(2)+Lfnorm(zind).*Jv(3); 

    f1=figure;
    xvals=src.pos(:,1);
    plot(xvals,Lfnorm_ori,'o','MarkerFaceColor',[0 .3 .3],'MarkerEdgeColor',[0 .5 .5],'MarkerSize',5);
    ax = gca;
    ax.FontSize = 18;
    box off
    set(gcf, 'Position', [570 769 670 227]);
    xlabel('Cranial-caudal position (cm)')
    ylabel('Leadfield norm')

    % Calculate the median y values for each unique x value
    %unique_x = unique(xvals);
   % median_y_values = splitapply(@median, Lfnorm_ori, findgroups(xvals));

    %hold on; % Hold the current plot
    %plot(unique_x(1:end-1), median_y_values(1:end-1), 'color',[0 .5 .5], 'LineWidth', 2);

% 
%     figure;
%     plot_func_spine_dat(subject,src,Lfnorm_ori,grad,sens_stl_ft)
%     title('leadfield strength')
    