function [coh_complex, wa, wb, ta, tb] = canonical_coherence(Cx, Cy, Cxy)
    % Cx: covariance matrix of spinal cord signals
    % Cy: covariance matrix of EMG signals
    % Cxy: cross-covariance matrix between spinal cord and EMG signals

    % Compute the inverse of Cy and square root of inverse of Cx
    Cy_inv = inv(Cy);
    Cx_invsqrt = sqrtm(inv(Cx));

    % Number of iterations for phase angle optimization
    n_iter = 10;
    n_iter1 = 5;
    coh = zeros(n_iter1, 1);
    
    % Find maximum coherence over different phase angles
    for n_iter = 1:n_iter1
        phi = n_iter / n_iter1 * pi;  % Phase angle
        coh(n_iter) = cohmax(phi, Cxy, Cy_inv, Cx_invsqrt);
    end
    
    % Find the phase angle that maximizes coherence
    [~, ind_cohmax] = max(coh);
    phi = ind_cohmax / n_iter1 * pi;
    
    % Newton's method to refine the optimal phase
    dphi = 0.000001;
    mu = 1;
    coh_0 = coh(ind_cohmax);
    for n_iter = 1:n_iter
        coh_0plus = cohmax(phi + dphi, Cxy, Cy_inv, Cx_invsqrt);
        coh_0minus = cohmax(phi - dphi, Cxy, Cy_inv, Cx_invsqrt);
        fprime = (coh_0plus - coh_0minus) / (2 * dphi);
        fprimeprime = (coh_0plus + coh_0minus - 2 * coh_0) / (dphi^2);
        deltaphi = -fprime / (fprimeprime - mu);
        phi_new = phi + deltaphi;
        phi_new = mod(phi_new + pi / 2, pi) - pi / 2;  % Ensure phase angle is within [-pi/2, pi/2]
        coh_0_new = cohmax(phi_new, Cxy, Cy_inv, Cx_invsqrt);
        
        if coh_0_new > coh_0
            mu = mu / 2;
            phi = phi_new;
            coh_0 = coh_0_new;
        else
            mu = mu * 2;
        end
    end
    
    % Compute the spatial weights that maximize coherence
    [coh_real, wa, wb] = cohmax_withdirs(phi, Cxy, Cy_inv, Cx_invsqrt);

    % Compute the canonical coherence for the complex values
    coh_complex = (wa' * Cxy * wb) / sqrt((wa' * Cx * wa) * (wb' * Cy * wb));
    
    % Compute the signal projections onto canonical variates
    ta = real(Cx) * real(wa);
    tb = real(Cy) * real(wb);
end

% Sub-function: _cohmax
function coh_val = cohmax(phi, Cxy, Cy_inv, Cx_invsqrt)
    % Compute the maximum coherence at a given phase angle phi
    Cxy = real(exp(-1i * phi) * Cxy);  % Apply the phase shift to Cxy
    X = Cxy * Cy_inv * Cxy';  % Cross-spectrum matrix
    Y = Cx_invsqrt * X * Cx_invsqrt';  % Transform into the Cx domain
    [~, S, ~] = svd(Y);  % Singular value decomposition
    coh_val = sqrt(S(1));  % The first singular value represents the coherence
end

% Sub-function: _cohmax_withdirs
function [coh_real, wa, wb] = cohmax_withdirs(phi, Cxy, Cy_inv, Cx_invsqrt)
    % Compute coherence and the spatial weights wa, wb at the given phase angle
    Cxy = real(exp(-1i * phi) * Cxy);  % Apply the phase shift to Cxy
    Cy_invsqrt = sqrtm(Cy_inv);  % Compute the square root of Cy_inv
    D = Cx_invsqrt * Cxy * Cy_invsqrt;  % Compute the intermediate matrix D
    Y = D * D';  % Compute Y matrix
    [U, S, ~] = svd(Y);  % Singular value decomposition
    wap = U(:, 1);  % First principal component of Y
    
    % Compute the spatial weights for both spinal cord and EMG signals
    wa = Cx_invsqrt * wap;  % Weight for spinal cord
    wa = wa / norm(wa);  % Normalize the weight
    
    wbp = D' * wap;  % EMG component corresponding to wap
    wb = Cy_invsqrt * wbp;  % Weight for EMG signal
    wb = wb / norm(wb);  % Normalize the weight
    
    coh_real = sqrt(S(1));  % The coherence value is the first singular value
end
