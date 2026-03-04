function plot_vip_bars(vip, pls_results, varargin)
% PLOT_VIP_BARS  Publication-quality VIP bar plots, colored by sensor.
%
%   Plots VIP>1 features as horizontal bars, colored by sensor identity,
%   with AEN and NEJ side by side. Includes VIP=1 threshold and model
%   summary in title.
%
% INPUTS
%   vip          - Struct with fields .AEN and/or .NEJ, each containing:
%                    .sorted  [m x 1] - sorted VIP scores (descending)
%                    .parms   {m x 1} - corresponding feature names
%   pls_results  - Output of fit_pls_model (for R², ncomp in title)
%
% OPTIONAL (name-value via ct_input)
%   cell_lines   - Cell lines to plot (default: fieldnames of vip)
%   max_features - Max features to show per panel (default: all VIP>1)

% --- Parse options ---
op.cell_lines   = fieldnames(vip)';
op.max_features = [];
op = ct_input(varargin, op);

% --- Sensor color map ---
cmap = struct( ...
    'ampk', [0.90 0.50 0.10], ... % orange
    'erkt', [0.85 0.33 0.10], ... % red-orange
    'ekar', [0.20 0.63 0.17], ... % green
    'jnkt', [0.58 0.40 0.74], ... % purple
    'nfkb', [0.22 0.47 0.78], ... % blue
    'nnuc', [0.30 0.68 0.87], ... % light blue
    'nnfc', [0.12 0.32 0.55]);    % dark blue

cls = op.cell_lines;
ncl = numel(cls);

figure('Position', [100 100 600*ncl 450], 'Name', 'VIP Scores');

for ci = 1:ncl
    cl = cls{ci};
    v  = vip.(cl);

    % Limit features if requested
    nf = numel(v.sorted);
    if ~isempty(op.max_features); nf = min(nf, op.max_features); end
    scores = v.sorted(nf:-1:1);     % Flip for bottom-to-top plotting
    parms  = v.parms(nf:-1:1);

    % Extract sensor tag per feature
    tags = regexp(parms, '_([a-z]+)$', 'tokens');
    tags = cellfun(@(x) x{1}{1}, tags, 'UniformOutput', false);

    % Assign colors
    colors = zeros(nf, 3);
    for fi = 1:nf
        if isfield(cmap, tags{fi})
            colors(fi,:) = cmap.(tags{fi});
        else
            colors(fi,:) = [0.5 0.5 0.5]; % fallback gray
        end
    end

    % Plot
    subplot(1, ncl, ci); hold on;
    for fi = 1:nf
        barh(fi, scores(fi), 'FaceColor', colors(fi,:), 'EdgeColor', 'none');
    end
    xline(1, '--k', 'LineWidth', 1);  % VIP=1 threshold

    yticks(1:nf);
    yticklabels(strrep(parms, '_', ' '));
    set(gca, 'FontSize', 8, 'TickLabelInterpreter', 'none');
    xlabel('VIP Score');
    ylim([0.5 nf+0.5]);

    % Title with model summary
    r2 = pls_results.(cl).R2;
    nc = pls_results.(cl).ncomp;
    title(sprintf('%s — VIP > 1  (ncomp=%d, R²=%.3f)', cl, nc, r2), ...
          'FontSize', 11);
    hold off;
end

% --- Build legend from all sensors across all panels ---
all_sensor_tags = {};
for ci = 1:ncl
    cl = cls{ci};
    v  = vip.(cl);
    nf = numel(v.sorted);
    if ~isempty(op.max_features); nf = min(nf, op.max_features); end
    parms = v.parms(1:nf);
    tags  = regexp(parms, '_([a-z]+)$', 'tokens');
    all_sensor_tags = [all_sensor_tags, cellfun(@(x) x{1}{1}, tags, 'UniformOutput', false)];
end
utags = unique(all_sensor_tags, 'stable');
h = gobjects(numel(utags), 1);
for ti = 1:numel(utags)
    h(ti) = patch(NaN, NaN, cmap.(utags{ti}), 'EdgeColor', 'none');
end
legend(h, utags, 'Location', 'southeast', 'FontSize', 8);

end