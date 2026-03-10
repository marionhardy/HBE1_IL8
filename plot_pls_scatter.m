function plot_pls_scatter(pls_results, varargin)
% PLOT_PLS_SCATTER  Predicted vs measured IL-8 from PLS model (Figure 3E).
%
%   Density scatter of PLS-predicted IL-8 vs measured IL-8 with
%   perfect prediction line (red dashed).
%
% INPUTS
%   pls_results  - Output of fit_pls_model
%
% OPTIONAL (name-value via ct_input)
%   cell_lines   - Cell lines to plot (default: fieldnames of pls_results)

op.cell_lines = fieldnames(pls_results)';
op = ct_input(varargin, op);

cls = op.cell_lines;
ncl = numel(cls);

figure('Position', [100 100 450*ncl 400], 'Name', 'Predicted vs Measured IL-8', 'Color', 'w');

for ci = 1:ncl
    cl = cls{ci};
    z = pls_results.(cl).pls_out;

    % Predicted Y = X * BETA(2:end) + BETA(1)
    y_pred = z.X * z.BETA(2:end) + z.BETA(1);
    y_meas = z.Y;

    subplot(1, ncl, ci);
    dscatter(y_meas, y_pred);
    hold on;

    % Perfect prediction line
    lims = [min([y_meas; y_pred]), max([y_meas; y_pred])];
    plot(lims, lims, 'r--', 'LineWidth', 1.5);

    % R²
    R2 = pls_results.(cl).R2;

    xlabel('Measured IL-8', 'FontSize', 10);
    ylabel('Predicted IL-8', 'FontSize', 10);
    title(sprintf('%s Cell Line  (R²=%.3f)', cl, R2), 'FontSize', 11);
    xlim(lims); ylim(lims);
    axis square; box on;
    hold off;
end
end