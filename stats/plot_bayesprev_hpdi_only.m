function out = plot_bayesprev_hpdi_only(indsig, alpha, titleStr)
% Minimal prevalence plot: MAP dot + 95% HPDI bar (no posterior curve)

if nargin < 2 || isempty(alpha)
    alpha = 0.05;
end
if nargin < 3
    titleStr = '';
end

indsig = logical(indsig(:));
k = sum(indsig);
n = numel(indsig);

xmap   = bayesprev_map(k,n,alpha);
hpdi95 = bayesprev_hpdi(0.95,k,n,alpha);

figure('Color','w','Position',[100 100 380 140]); hold on;

% 95% HPDI as a thick horizontal line
plot(hpdi95, [0 0], 'k-', 'LineWidth', 6);

% MAP as a dot
plot(xmap, 0, 'ko', 'MarkerFaceColor','k', 'MarkerSize',7);

xlim([0 1]);
ylim([-1 1]);
set(gca,'YTick',[],'TickDir','out','FontSize',12,'LineWidth',1.2);
xlabel('Population prevalence');

if ~isempty(titleStr)
    title(titleStr);
end
box off;

% Optional text inside (compact)
txt = sprintf('k=%d/%d, MAP=%.2f, 95%% HPDI=[%.2f, %.2f]', ...
    k, n, xmap, hpdi95(1), hpdi95(2));
text(0.01, 0.85, txt, 'Units','normalized', 'FontSize',12);

out.k = k; out.n = n; out.alpha = alpha; out.map = xmap; out.hpdi95 = hpdi95;
end