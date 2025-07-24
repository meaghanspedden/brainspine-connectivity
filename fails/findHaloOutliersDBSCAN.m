


  function outlierIdx = findHaloOutliersDBSCAN(haloPos, epsilon, minPoints)
    % haloPos: n x 3 matrix of positions
    % epsilon: the maximum distance between two points to be considered neighbors
    % minPoints: the minimum number of points required to form a dense region
    haloPos=[positions.Px positions.Py positions.Pz];
    minPoints=4;
    epsilon=120;

    % Apply DBSCAN to find outliers
    idx = dbscan(haloPos, epsilon, minPoints);
    
    % Identify the points labeled as noise (outliers are labeled as -1)
    outlierIdx = find(idx == -1);
    
    % Plot the result
    figure; hold on; grid on;
    scatter3(haloPos(:,1), haloPos(:,2), haloPos(:,3), 36, idx, 'filled'); % Clustered points
    scatter3(haloPos(outlierIdx,1), haloPos(outlierIdx,2), haloPos(outlierIdx,3), 'r', 'filled'); % Outliers
    xlabel('X'); ylabel('Y'); zlabel('Z');
    title('Outliers Detected Using DBSCAN');
    legend('Clusters', 'Outliers');
    
    % Display outlier indices
    disp('Outlier Indices:');
    disp(outlierIdx);

