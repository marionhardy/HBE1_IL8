function plot_sensor_scatter(S, sig_means, varargin)
% PLOT_SENSOR_SCATTER  Single sensor vs IL-8 density scatter (Figure 4B).
%
%   Density scatter of per-cell mean signal activity vs IL-8 intensity
%   with fitted linear model and R².
%
% INPUTS
%   S           - Main analysis struct
%   sig_means   - Output of extract_signal_means_from_d
%
% OPTIONAL (name-value via ct_input)
%   signals_AEN - Signals to plot for AEN
%   signals_NEJ - Signals to plot for NEJ
%   cell_lines  - Cell lines to plot (default: {'AEN','NEJ'})

op.signals_aen = {'Mean_erktr','Mean_ampkar','Mean_nfkb_n'};
op.signals_nej = {'Mean_ekar','Mean_jnktr','Mean_nfkb_n'};
op.cell_lines  = {'AEN','NEJ'};
op = ct_input(varargin, op);

sig_field = struct('AEN', {op.signals_aen}, 'NEJ', {op.signals_nej});
cls = op.cell_lines;

for ci = 1:numel(cls)
    cl = cls{ci};
    sigs = sig_field.(cl);
    nsig = numel(sigs);
    idx = S.idx.(cl);
    y = S.Y.IL8c_raw(idx);

    figure('Position', [100 100 350*nsig 350], ...
           'Name', sprintf('%s Sensor vs IL-8', cl), 'Color', 'w');

    for si = 1:nsig
        x = sig_means{idx, sigs{si}};

        % Remove NaN
        valid = ~isnan(x) & ~isnan(y);
        xv = x(valid);
        yv = y(valid);

        subplot(1, nsig, si);
        % Add small jitter to prevent dscatter banding
        jit = 0.05 * std(yv);
        dscatter(xv, yv + jit*randn(size(yv)));
        hold on;

        % Linear fit
        p = polyfit(xv, yv, 1);
        xlims = [min(xv), max(xv)];
        xfit = linspace(xlims(1), xlims(2), 100)';
        yfit = polyval(p, xfit);
        plot(xfit, yfit, 'r-', 'LineWidth', 1.5);

        % R²
        SS_res = sum((yv - polyval(p, xv)).^2);
        SS_tot = sum((yv - mean(yv)).^2);
        R2 = 1 - SS_res / SS_tot;

        % Axis limits: full lower range, clip upper at 99.75th percentile
        pct = 99.75;
        xlim([min(xv), prctile(xv, pct)]);
        ylim([min(yv), prctile(yv, pct)]);

        % Labels
        sig_label = strrep(strrep(sigs{si}, 'Mean_', ''), '_', ' ');
        xlabel(sig_label, 'FontSize', 10);
        ylabel('IL-8 Intensity', 'FontSize', 10);
        title(sprintf('R²: %.3f', R2), 'FontSize', 11);
        axis square; box on;
        hold off;
    end

    sgtitle(sprintf('%s Cell Line', cl), 'FontSize', 12);
end
end