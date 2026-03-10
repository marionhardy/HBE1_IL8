function [results] = fit_reduced_pls(S, varargin)
% FIT_REDUCED_PLS  Bootstrap VIP selection + CV-evaluated reduced model.
%
%   Step 1: Bootstrap VIP stability — resample cells with replacement,
%           fit PLS, compute VIP scores. Keep features with VIP>1 in
%           >= stability_thresh fraction of bootstraps.
%   Step 2: Fit reduced model (stable features only) with k-fold CV.
%   Step 3: Compare full vs reduced model performance.
%
% INPUTS
%   S               - Main analysis struct with S.pls.(cl).X, .parms, .y
%
% OPTIONAL (name-value via ct_input)
%   ncomp           - Components for full model (default: 10)
%   cv              - CV folds for final evaluation (default: 10)
%   n_boot          - Number of bootstrap iterations (default: 500)
%   stability_thresh - Fraction of bootstraps feature must have VIP>1 (default: 0.50)
%   cell_lines      - Cell lines to run (default: {'AEN','NEJ'})
%   rm_xsout        - Pass to pls.m (default: false)
%   plot            - Generate diagnostic plots (default: true)
%
% OUTPUTS
%   results         - Struct per cell line:
%       .stable_parms   - Features passing stability threshold
%       .vip_freq       [1 x nFeatures] - Fraction of bootstraps with VIP>1
%       .full_R2        - Full model R²
%       .full_cv_mse    - Full model CV MSE
%       .reduced_R2     - Reduced model R²
%       .reduced_cv_mse - Reduced model CV MSE
%       .full_pls       - Full model pls.m output
%       .reduced_pls    - Reduced model pls.m output
%       .n_cells_full   - Cells retained in full model
%       .n_cells_reduced - Cells retained in reduced model

% --- Parse options ---
op.ncomp            = 10;
op.cv               = 10;
op.n_boot           = 500;
op.stability_thresh = 0.50;
op.cell_lines       = {'AEN','NEJ'};
op.rm_xsout         = false;
op.plot             = true;
op = ct_input(varargin, op);

cls = op.cell_lines;
results = struct();

for ci = 1:numel(cls)
    cl = cls{ci};
    fprintf('\n=== %s: Bootstrap VIP stability (%d iterations) ===\n', cl, op.n_boot);

    % --- Extract data ---
    parms = S.pls.(cl).parms;
    X     = S.pls.(cl).X{:, parms};
    y     = S.pls.(cl).y;
    [N, nFeat] = size(X);
    fprintf('  Cells: %d | Features: %d\n', N, nFeat);

    % --- Bootstrap VIP ---
    vip_count = zeros(1, nFeat);  % Count of VIP>1 per feature

    for bi = 1:op.n_boot
        % Resample with replacement
        idx = randi(N, N, 1);
        X_boot = X(idx, :);
        y_boot = y(idx);

        % Fit PLS
        z_boot = pls(X_boot, y_boot, 'ncomp', op.ncomp, 'cv', 'resubstitution', ...
                     'ploton', false, 'params', parms, 'append', true, ...
                     'rm_xsout', op.rm_xsout);

        % Compute VIP
        [~, ~, vip_all] = pls_VIP_score(z_boot, parms);
        vip_count = vip_count + (vip_all' > 1);

        if mod(bi, 100) == 0
            fprintf('  Bootstrap %d/%d\n', bi, op.n_boot);
        end
    end

    % --- Stability selection ---
    vip_freq = vip_count / op.n_boot;
    stable_mask = vip_freq >= op.stability_thresh;
    stable_parms = parms(stable_mask);

    fprintf('  Stable features (VIP>1 in >= %.0f%% of bootstraps): %d / %d\n', ...
        op.stability_thresh * 100, sum(stable_mask), nFeat);
    fprintf('  '); disp(stable_parms');

    % --- Fit full model with CV ---
    fprintf('  Fitting full model (CV=%d)...\n', op.cv);
    z_full = pls(X, y, 'ncomp', op.ncomp, 'cv', op.cv, 'ploton', false, ...
                 'params', parms, 'append', true, 'rm_xsout', op.rm_xsout);
    full_R2     = sum(z_full.PCTVAR(2,:));
    full_cv_mse = z_full.MSE(2, end);

    % --- Fit reduced model with CV ---
    nc_reduced = min(op.ncomp, sum(stable_mask));
    fprintf('  Fitting reduced model (%d features, ncomp=%d, CV=%d)...\n', ...
        sum(stable_mask), nc_reduced, op.cv);
    z_reduced = pls(X(:, stable_mask), y, 'ncomp', nc_reduced, 'cv', op.cv, ...
                    'ploton', false, 'params', stable_parms, 'append', true, ...
                    'rm_xsout', op.rm_xsout);
    reduced_R2     = sum(z_reduced.PCTVAR(2,:));
    reduced_cv_mse = z_reduced.MSE(2, end);

    % --- Report ---
    fprintf('\n  %s Results:\n', cl);
    fprintf('    Full model:    R²=%.4f  CV MSE=%.4f  n=%d  (%d features)\n', ...
        full_R2, full_cv_mse, size(z_full.X,1), nFeat);
    fprintf('    Reduced model: R²=%.4f  CV MSE=%.4f  n=%d  (%d features)\n', ...
        reduced_R2, reduced_cv_mse, size(z_reduced.X,1), sum(stable_mask));

    % --- Store ---
    results.(cl).stable_parms    = stable_parms;
    results.(cl).stable_mask     = stable_mask;
    results.(cl).vip_freq        = vip_freq;
    results.(cl).full_R2         = full_R2;
    results.(cl).full_cv_mse     = full_cv_mse;
    results.(cl).reduced_R2      = reduced_R2;
    results.(cl).reduced_cv_mse  = reduced_cv_mse;
    results.(cl).full_pls        = z_full;
    results.(cl).reduced_pls     = z_reduced;
    results.(cl).n_cells_full    = size(z_full.X, 1);
    results.(cl).n_cells_reduced = size(z_reduced.X, 1);
    results.(cl).parms           = parms;
end

% --- Plot ---
if op.plot
    plot_vip_stability(results, cls, op.stability_thresh);
end

end


%% ========================================================================
function plot_vip_stability(results, cls, thresh)
% Bootstrap VIP frequency per feature, colored by sensor.

cmap = struct( ...
    'ampk', [0.90 0.50 0.10], ...
    'erkt', [0.85 0.33 0.10], ...
    'ekar', [0.20 0.63 0.17], ...
    'jnkt', [0.58 0.40 0.74], ...
    'nfkb', [0.22 0.47 0.78], ...
    'nnuc', [0.30 0.68 0.87], ...
    'nnfc', [0.12 0.32 0.55]);

ncl = numel(cls);
figure('Position', [100 100 600*ncl 500], 'Name', 'VIP Stability');

for ci = 1:ncl
    cl = cls{ci};
    r  = results.(cl);

    % Sort by frequency
    [freq_sorted, si] = sort(r.vip_freq, 'descend');
    parms_sorted = r.parms(si);

    % Only show features that ever crossed VIP>1
    show = freq_sorted > 0;
    freq_sorted = freq_sorted(show);
    parms_sorted = parms_sorted(show);
    nf = numel(freq_sorted);

    % Sensor tags + colors
    tags = regexp(parms_sorted, '_([a-z]+)$', 'tokens');
    tags = cellfun(@(x) x{1}{1}, tags, 'UniformOutput', false);
    colors = zeros(nf, 3);
    for fi = 1:nf
        if isfield(cmap, tags{fi})
            colors(fi,:) = cmap.(tags{fi});
        else
            colors(fi,:) = [0.5 0.5 0.5];
        end
    end

    % Plot (bottom to top)
    subplot(1, ncl, ci); hold on;
    for fi = 1:nf
        barh(nf - fi + 1, freq_sorted(fi), 'FaceColor', colors(fi,:), 'EdgeColor', 'none');
    end
    xline(thresh, '--k', 'LineWidth', 1.2);
    yticks(1:nf);
    yticklabels(strrep(parms_sorted(nf:-1:1), '_', ' '));
    set(gca, 'FontSize', 7, 'TickLabelInterpreter', 'none');
    xlabel('Fraction of bootstraps with VIP > 1');
    title(sprintf('%s — VIP Stability (full: R²=%.3f, reduced: R²=%.3f)', ...
        cl, r.full_R2, r.reduced_R2), 'FontSize', 10);
    xlim([0 1]);
    hold off;

    % Legend
    utags = unique(tags, 'stable');
    h = gobjects(numel(utags), 1);
    for ti = 1:numel(utags)
        h(ti) = patch(NaN, NaN, cmap.(utags{ti}), 'EdgeColor', 'none');
    end
    legend(h, utags, 'Location', 'southeast', 'FontSize', 7);
end

sgtitle('Bootstrap VIP Feature Stability');
end