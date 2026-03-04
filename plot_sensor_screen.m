function plot_sensor_screen(screen, varargin)
% PLOT_SENSOR_SCREEN  Visualize sensor combination screen results.
%
%   Plots CV MSE by sensor combination, colored by number of sensors,
%   with AEN and NEJ side by side.
%
% INPUTS
%   screen       - Output of run_sensor_screen (struct with .AEN, .NEJ tables)
%
% OPTIONAL (name-value via ct_input)
%   cell_lines   - Cell lines to plot (default: fieldnames of screen)
%   sort_by      - 'cv_mse' or 'R2' (default: 'cv_mse')
%   show_perm    - Mark models with significant permutation p (default: true)

op.cell_lines = fieldnames(screen)';
op.sort_by    = 'cv_mse';
op.show_perm  = true;
op = ct_input(varargin, op);

cls = op.cell_lines;
ncl = numel(cls);
figure('Position', [100 100 600*ncl 500], 'Name', 'Sensor Screen');

for ci = 1:ncl
    cl = cls{ci};
    T  = screen.(cl);

    % Sort by number of sensors, then by sort metric within group
    n_sens = cellfun(@(x) numel(strsplit(x, '+')), T.sensors);
    [~, si] = sortrows([n_sens, T.(op.sort_by)], [1, 2]);
    T = T(si, :);
    n_sens = n_sens(si);

    max_sens = max(n_sens);
    cmap = parula(max_sens + 1);

    subplot(1, ncl, ci); hold on;
    for ri = 1:height(T)
        barh(ri, T.(op.sort_by)(ri), 'FaceColor', cmap(n_sens(ri),:), 'EdgeColor', 'none');
    end

    % Mark significant permutation results
    if op.show_perm && ismember('perm_p', T.Properties.VariableNames)
        sig_idx = find(T.perm_p <= 0.05 & ~isnan(T.perm_p));
        for si2 = 1:numel(sig_idx)
            ri = sig_idx(si2);
            text(T.(op.sort_by)(ri) + 0.005, ri, '*', ...
                 'FontSize', 14, 'FontWeight', 'bold', 'Color', 'r');
        end
    end

    yticks(1:height(T));
    yticklabels(strrep(T.sensors, '_', ' '));
    set(gca, 'FontSize', 7, 'TickLabelInterpreter', 'none');
    xlabel(strrep(op.sort_by, '_', ' '));
    title(sprintf('%s — Sensor Combinations', cl), 'FontSize', 11);
    ylim([0.5 height(T)+0.5]);

    % Legend for sensor count
    h = gobjects(max_sens, 1);
    for k = 1:max_sens
        h(k) = patch(NaN, NaN, cmap(k,:), 'EdgeColor', 'none');
    end
    labels = arrayfun(@(x) sprintf('%d sensor(s)', x), 1:max_sens, ...
                      'UniformOutput', false);
    if op.show_perm
        h(end+1) = plot(NaN, NaN, 'r*', 'MarkerSize', 10);
        labels{end+1} = 'perm p \leq 0.05';
    end
    legend(h, labels, 'Location', 'southeast', 'FontSize', 7);
    hold off;
end

sgtitle('Sensor Combination Screen');
end