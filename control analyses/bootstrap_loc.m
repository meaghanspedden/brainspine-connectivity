clear all; close all; clc

subs = {'OP00212', 'OP00213',  'OP00215', 'OP00219', ...
        'OP00220', 'OP00221', 'OP00224', 'OP00225', 'OP00226'};
which_ori = 'all';  % 
HFC = 1;
gemfile = 'D:\MSST001\generic_merged\geoms.mat';
castforward = load(gemfile);
src=castforward.sources_center_line;


cmap = lines(length(subs));  
y_meds=[];
mad_y=[];

for i = 1:length(subs)
    sub = subs{i};
    save_dir = fullfile('D:\MSST001', [sub '_contrast']);
    
    if HFC
        savename = sprintf('spine_bootstrap_%s_%s.mat', sub, which_ori);
    else
        savename = sprintf('spine_bootstrap_%s_%s_nohfc.mat', sub, which_ori);
    end

    data = load(fullfile(save_dir, savename));
    y_pos = data.bootstrap_y_positions;

    y_meds(i)=data.med_y;
    mad_y(i) = mad(y_pos, 1);


end

% xlabel('Spinal cord Y-position (mm)');
% 
% legend('show');

%%
num_subjects = length(subs);
colors = lines(num_subjects);
x_offset = .2;  % Spacing between subjects
ycoord = src.pos(:,2);

figure; hold on;

for subj = 1:num_subjects
    % Convert MAD to approximate standard deviation
    mad_std = mad_y(subj) * 1.4826;
    if mad_std == 0
    mad_std = 1e-6;  % Small fallback value
    end
    med = y_meds(subj);

    % Gaussian PDF
    gauss_pdf = (1 / (mad_std * sqrt(2*pi))) * exp(-0.5 * ((ycoord - med)/mad_std).^2);

    % Scale if desired (e.g., for visual consistency)
    gauss_pdf_scaled = gauss_pdf / max(gauss_pdf);  % optional: normalize for display
    gauss_pdf_scaled = gauss_pdf_scaled * 0.4;       % control width for visibility

    x_pos = subj * x_offset;

    % Plot sideways Gaussian
    plot(gauss_pdf_scaled + x_pos, ycoord, 'Color', colors(subj,:), 'LineWidth', 2);

    % Median line
%     plot([x_pos - 0.2, x_pos + 0.2], [med, med], 'k-', 'LineWidth', 2);
% 
%     % ±2 MAD shading
%     lower = med - 2 * mad_y(subj);
%     upper = med + 2 * mad_y(subj);
%     fill([x_pos-0.1, x_pos+0.1, x_pos+0.1, x_pos-0.1], ...
%          [lower, lower, upper, upper], ...
%          colors(subj,:), 'FaceAlpha', 0.2, 'EdgeColor', 'none');
%     disp('done')
end

xlabel('Subject');
ylabel('Spinal cord Y-position (mm)');
title('Median ± 2 MAD');
xlim([0, (num_subjects+1)*x_offset]);
set(gca, 'XTick', x_offset*(1:num_subjects), 'XTickLabel', 1:num_subjects);

%%
colors = parula(num_subjects);
ycoord = src.pos(:,2);         % Y-axis coordinate
x_pos = 1;                     % Shared X-position for all

figure; hold on;

for subj = 1:num_subjects
    med = y_meds(subj);
    mad_std = mad_y(subj) * 1.4826;

    % Handle degenerate (zero spread) case
    if mad_std == 0
        mad_std = 1e-6;
    end

    % Gaussian PDF over ycoord
    gauss_pdf = (1 / (mad_std * sqrt(2*pi))) * exp(-0.5 * ((ycoord - med)/mad_std).^2);

    % Normalize and scale for visual appearance
    gauss_pdf_scaled = gauss_pdf / max(gauss_pdf);  % Normalize
    gauss_pdf_scaled = gauss_pdf_scaled * 0.4;       % Control width

    % Plot Gaussian sideways at same x_pos
    plot(gauss_pdf_scaled + x_pos, ycoord, 'Color', colors(subj,:), 'LineWidth', 1.5);
    disp('done')
%     % Optional: plot median line
%     plot([x_pos - 0.1, x_pos + 0.1], [med, med], 'k--', 'LineWidth', 1);
% 
%     % Optional: ±2 MAD shading
%     lower = med - 2 * mad_y(subj);
%     upper = med + 2 * mad_y(subj);
%     fill([x_pos - 0.1, x_pos + 0.1, x_pos + 0.1, x_pos - 0.1], ...
%          [lower, lower, upper, upper], ...
%          colors(subj,:), 'FaceAlpha', 0.15, 'EdgeColor', 'none');
end

% Final plot formatting
xlabel('Distribution');
ylabel('Spinal cord y-position (mm)');
title('Median ± 2 MAD');
xlim([x_pos - 0.2, x_pos + 0.6]);
set(gca, 'XTick', []);
box on;
legend(subs)