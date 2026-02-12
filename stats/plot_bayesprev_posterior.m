function out = plot_bayesprev_posterior(indsig, alpha, titleStr)
% Publication-ready Bayesian population prevalence plot
% with informative annotation (k/n, alpha, MAP)
%
% indsig  : logical or 0/1 vector (length n)
% alpha   : false positive rate (default 0.05)
% titleStr: optional title ('' for none)

if nargin < 2 || isempty(alpha)
    alpha = 0.05;
end
if nargin < 3
    titleStr = '';
end

indsig = logical(indsig(:));
k = sum(indsig);
n = numel(indsig);

% Grid for posterior
x = linspace(0,1,400);
posterior = bayesprev_posterior(x,k,n,alpha);

% MAP
xmap = bayesprev_map(k,n,alpha);
pmap = bayesprev_posterior(xmap,k,n,alpha);

% Credible intervals
hpdi95 = bayesprev_hpdi(0.95,k,n,alpha);
hpdi50 = bayesprev_hpdi(0.50,k,n,alpha);

% Lower bound
bound = bayesprev_bound(0.95,k,n,alpha);

%% ---- Plot ----
figure('Color','w','Position',[100 100 600 450]); hold on;

% Posterior curve
plot(x, posterior, 'k', 'LineWidth', 2);

% 95% HPDI (light grey)
plot(hpdi95, [pmap pmap], ...
    'Color',[0.6 0.6 0.6], ...
    'LineWidth',2);

% 50% HPDI (darker grey)
plot(hpdi50, [pmap pmap], ...
    'Color',[0.3 0.3 0.3], ...
    'LineWidth',4);

% MAP point
plot(xmap, pmap, 'ko', ...
    'MarkerFaceColor','k', ...
    'MarkerSize',7);

% Lower 95% bound
line([bound bound], [0 bayesprev_posterior(bound,k,n,alpha)], ...
    'Color',[0.4 0.4 0.4], ...
    'LineStyle','--', ...
    'LineWidth',1.5);

% Labels
xlabel('Population prevalence proportion');
ylabel('Posterior density');

if ~isempty(titleStr)
    title(titleStr);
end

set(gca, ...
    'FontSize',14, ...
    'LineWidth',1.2, ...
    'TickDir','out');

box off;

% Informative annotation (top-left inside panel)
annotationText = sprintf('k = %d / %d\n\\alpha = %.3f\nMAP = %.2f\n95%% HPDI = [%.2f, %.2f]', ...
    k, n, alpha, xmap, hpdi95(1), hpdi95(2));

text(0.02, 0.98, annotationText, ...
    'Units','normalized', ...
    'VerticalAlignment','top', ...
    'FontSize',14);

%% ---- Output structure ----
out.k = k;
out.n = n;
out.alpha = alpha;
out.map = xmap;
out.bound95 = bound;
out.hpdi95 = hpdi95;
out.hpdi50 = hpdi50;
end