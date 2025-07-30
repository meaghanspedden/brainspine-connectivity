function smoothed_outline = smooth_head_outline_spline(XZ, npoints, smooth_factor)
    % Inputs:
    %   XZ           - [N x 2] head outline in X-Z coordinates
    %   npoints      - number of output points (default: 200)
    %   smooth_factor - 0 (no smoothing) to 1 (max smoothing), default: 0.9

    if nargin < 2
        npoints = 200;
    end
    if nargin < 3
        smooth_factor = 0.9;
    end

    % Ensure closed loop
    if ~isequal(XZ(1,:), XZ(end,:))
        XZ = [XZ; XZ(1,:)];  % close shape if not already
    end

    % Parametrize by arc length
    d = [0; cumsum(sqrt(sum(diff(XZ).^2,2)))];  % arc length
    t = linspace(0, d(end), npoints);           % uniform arc length

    % Fit smoothing splines to X and Z
    px = csaps(d, XZ(:,1), 1 - smooth_factor);
    pz = csaps(d, XZ(:,2), 1 - smooth_factor);

    Xs = fnval(px, t);
    Zs = fnval(pz, t);

    smoothed_outline = [Xs(:), Zs(:)];
end
