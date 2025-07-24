function outlierIdx = findHaloOutliers(haloPos)
    % This function lets the user select outliers from a 3D scatter plot.
    % haloPos is expected to be an n x 3 matrix (positions of sensors in 3D space).

haloPos=[positions.Px positions.Py positions.Pz];




    % Plot the sensor positions
    figure; hold on; grid on;
    scatter3(haloPos(:,1), haloPos(:,2), haloPos(:,3), 'b', 'filled'); % All sensors
    xlabel('X'); ylabel('Y'); zlabel('Z');
    title('Select Outliers by Clicking');
    axis equal;
    
    % Initialize output
    outlierIdx = [];

    % Define view angles (e.g., top view, side view, and another perspective)
    viewAngles = [
        0, 90;  % View from the top
        90, 0;  % Side view
        45, 45; % Another perspective
    ];

    % Loop through each view and allow the user to select outliers
    for i = 1:size(viewAngles, 1)
        % Set the view to the current angles
        view(viewAngles(i, 1), viewAngles(i, 2));
        
        % Allow the user to click and select points
        [x_selected, y_selected] = getpts;  % Click on the plot to select points

        % Loop through the selected points to identify the closest sensor positions
        for j = 1:length(x_selected)
            % Calculate the distance between the selected point and all sensor positions
            dists = sqrt((haloPos(:,1) - x_selected(j)).^2 + (haloPos(:,2) - y_selected(j)).^2);
            
            % Find the closest sensor position to the selected point
            [~, idx] = min(dists);
            
            % Add the index of the closest sensor to the list of outliers
            outlierIdx = unique([outlierIdx; idx]);
        end
    end
    
    % Mark the outliers on the plot
    scatter3(haloPos(outlierIdx,1), haloPos(outlierIdx,2), haloPos(outlierIdx,3), 'r', 'filled');
    
    % Display the selected outlier indices
    disp('Selected Outlier Indices:');
    disp(outlierIdx);
    
    % Add a legend to the plot
    legend('Sensors', 'Selected Outliers');
end
