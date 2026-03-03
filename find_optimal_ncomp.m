function [results] = find_optimal_ncomp(S, varargin)
% FIND_OPTIMAL_NCOMP  CV-based component selection for PLS regression.
%
%   Sweeps ncomp = 1:max_ncomp, calling pls.m at each step with k-fold CV.
%   Reports optimal ncomp (min CV MSE) and parsimonious ncomp (1-SE rule).
%   Captures resubstitution R² alongside CV MSE for reference.
%
% INPUTS
%   S           - Main analysis struct with S.pls.(cl).X, .parms, .y
%
% OPTIONAL (name-value via ct_input)
%   max_ncomp   - Maximum components to test (default: 20)
%   cv          - Number of CV folds (default: 10)
%   cell_lines  - Cell lines to run (default: {'AEN','NEJ'})
%   plot        - Generate diagnostic plots (default: true)
%
% OUTPUTS
%   results     - Struct with fields per cell line:
%       .cv_mse     [1 x max_ncomp]  - CV MSE for Y at each ncomp
%       .cv_mse0    [scalar]          - Intercept-only CV MSE
%       .R2_cum     [1 x max_ncomp]  - Cumulative resubstitution R²
%       .ncomp_opt  [scalar]          - ncomp at min CV MSE
%       .ncomp_1se  [scalar]          - Parsimonious ncomp (1-SE rule)
%       .se         [1 x max_ncomp]  - SE of CV MSE (estimated)
%       .n_cells    [1 x max_ncomp]  - Cells retained per ncomp (after outlier removal)

% --- Parse options ---
op.max_ncomp  = 20;
op.cv         = 10;
op.cell_lines = {'AEN','NEJ'};
op.plot       = true;
op = ct_input(varargin, op);

cls = op.cell_lines;
results = struct();

for ci = 1:numel(cls)
    cl = cls{ci};
    fprintf('\n=== %s: sweeping ncomp 1:%d ===\n', cl, op.max_ncomp);

    % --- Extract data ---
    X = S.pls.(cl).X{:, S.pls.(cl).parms};  % Named column indexing
    y = S.pls.(cl).y;
    [N, nFeat] = size(X);
    max_nc = min(op.max_ncomp, nFeat);
    fprintf('  Cells: %d | Features: %d | Max ncomp: %d\n', N, nFeat, max_nc);

    % --- Preallocate ---
    cv_mse  = nan(1, max_nc);
    cv_mse0 = nan;
    R2_cum  = nan(1, max_nc);
    n_cells = nan(1, max_nc);

    % --- Sweep ---
    for nc = 1:max_nc
        z = pls(X, y, 'ncomp', nc, 'cv', op.cv, 'ploton', false, ...
                'params', S.pls.(cl).parms, 'append', true);

        % z.MSE: [2 x (nc+1)], row 2 = Y, col 1 = intercept-only
        cv_mse0 = z.MSE(2, 1);         % Same regardless of nc
        cv_mse(nc) = z.MSE(2, end);    % CV MSE at this ncomp

        % Cumulative resubstitution R²
        R2_cum(nc) = sum(z.PCTVAR(2, :));

        % Track how many cells survived outlier removal
        n_cells(nc) = size(z.X, 1);

        fprintf('  ncomp=%2d | CV MSE=%.4f | R²=%.4f | n=%d\n', ...
                nc, cv_mse(nc), R2_cum(nc), n_cells(nc));
    end

    % --- Optimal ncomp: minimum CV MSE ---
    [~, ncomp_opt] = min(cv_mse);

    fprintf('  >> Optimal ncomp: %d (CV MSE = %.4f)\n', ncomp_opt, cv_mse(ncomp_opt));

    % --- Store results ---
    results.(cl).cv_mse    = cv_mse;
    results.(cl).cv_mse0   = cv_mse0;
    results.(cl).R2_cum    = R2_cum;
    results.(cl).ncomp_opt = ncomp_opt;
    results.(cl).n_cells   = n_cells;
end

% --- Plot ---
if op.plot
    plot_ncomp_results(results, cls);
end

end


%% ========================================================================
function plot_ncomp_results(results, cls)
% Side-by-side CV MSE (with n_cells overlay) and R² vs ncomp.

figure('Position', [100 100 1200 500], 'Name', 'find_optimal_ncomp');

for ci = 1:numel(cls)
    cl = cls{ci};
    r  = results.(cl);
    nc = 1:numel(r.cv_mse);

    % --- CV MSE panel with n_cells overlay ---
    subplot(2, numel(cls), ci);
    yyaxis left; hold on;
    plot(nc, r.cv_mse, 'ko-', 'MarkerFaceColor', 'k', 'MarkerSize', 4);

    % Mark optimal ncomp
    plot(r.ncomp_opt, r.cv_mse(r.ncomp_opt), 'rv', 'MarkerSize', 10, ...
         'MarkerFaceColor', 'r');

    % Intercept-only reference
    yline(r.cv_mse0, '--', 'Color', [0.5 0.5 0.5]);
    ylabel('CV MSE (Y)');
    hold off;

    yyaxis right;
    plot(nc, r.n_cells, '-', 'Color', [0.2 0.6 0.2], 'LineWidth', 1.2);
    ylabel('n cells retained');

    xlabel('ncomp');
    title(sprintf('%s — CV MSE', cl));
    legend({'CV MSE', 'Min CVE', 'Intercept-only', 'n cells'}, ...
           'Location', 'best', 'FontSize', 7);

    % --- R² panel ---
    subplot(2, numel(cls), numel(cls) + ci); hold on;
    plot(nc, r.R2_cum, 'ko-', 'MarkerFaceColor', 'k', 'MarkerSize', 4);
    plot(r.ncomp_opt, r.R2_cum(r.ncomp_opt), 'rv', 'MarkerSize', 10, ...
         'MarkerFaceColor', 'r');
    xlabel('ncomp'); ylabel('Cumulative R² (resub)');
    title(sprintf('%s — R² (resubstitution)', cl));
    legend({'R²', 'Min CVE'}, 'Location', 'best', 'FontSize', 7);
    hold off;
end

sgtitle('Optimal ncomp Selection');
end