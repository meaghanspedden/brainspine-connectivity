% === Parameters ===
freqband = [10 35]; 
tapsmofrq = 2;

n_orient = 3;

% === Step 1: Fourier analysis per orientation ===
for o = 1:n_orient
    for k = 1:size(st.wdatastatic,2)
        alldat.trial{k}(:,:) = squeeze(st.wdatastatic(:,k,:,o));
    end

    for j = 1:numel(alldat.label)
        alldat.label{j} = sprintf('source%g', j);
    end

    cfg = [];
    cfg.output = 'fourier';  % <-- complex output
    cfg.method = 'mtmfft';
    cfg.foilim = freqband;
    cfg.tapsmofrq = tapsmofrq;
    cfg.keeptrials = 'yes';
    cfg.pad = 'nextpow2';

    fstat_fourier{o} = ft_freqanalysis(cfg, alldat);
end

% === Step 2: Stack into [3 x freq x src x trial] matrix ===
n_src    = length(fstat_fourier{1}.label);
n_freq   = length(fstat_fourier{1}.freq);
n_trials = size(st.wdatastatic,2);

J = zeros(3, n_freq, n_src, n_trials);

for o = 1:n_orient
    % fourierspctrm: trials x sources x freq
    J(o, :, :, :) = permute(fstat_fourier{o}.fourierspctrm(1:ntrials,:,:), [4, 3, 2, 1]);
end

% === Step 3: Estimate dominant orientation ===
dominant_orient = zeros(3, n_freq, n_src);
J_projected = zeros(n_freq, n_src, n_trials);

for s = 1:n_src
    for f = 1:n_freq
        Jvec = squeeze(J(:, f, s, :)); % [3 x trials], complex
        Jcov = cov(real(Jvec)');       % You can also try cov(abs(Jvec)') if you prefer

        % Get direction of max variance
        [V, D] = eig(Jcov);
        [~, idx] = max(diag(D));
        v_opt = V(:, idx);
        v_opt = v_opt / norm(v_opt);

        dominant_orient(:, f, s) = v_opt;

        % Project each trial
        for tr = 1:n_trials
            J_projected(f, s, tr) = v_opt' * Jvec(:, tr);
        end
    end
end



dominant_vecs = mean(dominant_orient, 2);  % [3 x 1 x n_src]
dominant_vecs = squeeze(dominant_vecs);    % [3 x n_src]

% Normalize vectors for consistent arrow size
vecs = dominant_vecs';
vecs = vecs ./ vecnorm(vecs, 2, 2);  % [n_src x 3]

vecs(end,:)=[];

% Plot each vector at the source position
% figure;
% hold on
% ft_plot_mesh(mtorso,'FaceColor', [0.8 0.8 1.0], 'facealpha', 0.2, 'EdgeColor','none',...
%     'AmbientStrength', 0.15);
% quiver3(src.pos(:,1), src.pos(:,2), src.pos(:,3), ...
%         vecs(:,1), vecs(:,2), vecs(:,3), ...
%         0.8, 'LineWidth', 1.5, 'Color', [0.1 0.5 0.9]);
% axis equal; grid on;
% xlabel('X'); ylabel('Y'); zlabel('Z');
% 
% plot3(src.pos(peakind2use,1), src.pos(peakind2use,2), src.pos(peakind2use,3), ...
%     'o', 'MarkerSize', 12, 'MarkerFaceColor', 	[1 0.4 0], 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
% 
% title('Dominant Orientation Vectors per Source');
% 
% view(90,11)


