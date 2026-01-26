

nPos = numel(Lf.leadfield);   % e.g. 44
lf_strength = nan(nsourcepoints,1);

for i = 1:nsourcepoints
    L = Lf.leadfield{i};      % [Nchan x 3]
    lf_strength(i) = norm(L, 'fro');
end


x=sources_cent.pos(:,2); figure; hold on;
plot(x,lf_strength, 'o-','LineWidth',2, 'MarkerFaceColor', [0.6 0.8 1], 'MarkerSize', 8)
xlabel('Cranial caudal position (mm)');
ylabel('Leadfield strength')
set(gca, 'FontSize', 13)
grid on;