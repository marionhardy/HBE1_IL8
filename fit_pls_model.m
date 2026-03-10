function [results] = fit_pls_model(S, varargin)
% FIT_PLS_MODEL  Fit PLSR model per cell line with permutation test.
%
%   Calls pls.m with explicit ncomp and k-fold CV. Runs scramble control
%   by default to establish significance via permutation p-value.
%
% INPUTS
%   S               - Main analysis struct with S.pls.(cl).X, .parms, .y
%
% OPTIONAL (name-value via ct_input)
%   ncomp_results   - Output of find_optimal_ncomp; extracts ncomp per cell line
%   ncomp           - Manual override: scalar applied to all cell lines
%   cv              - Number of CV folds (default: 10)
%   cell_lines      - Cell lines to run (default: {'AEN','NEJ'})
%   scramble        - Run permutation test (default: true)
%   n_perm          - Number of permutations (default: 100)
%   use_1se         - Use 1-SE ncomp instead of optimal (default: false)
%
% OUTPUTS
%   results         - Struct with fields per cell line:
%       .pls_out        - Full pls.m output struct
%       .ncomp          - ncomp used
%       .R2             - Cumulative resubstitution R² (Y)
%       .cv_mse         - CV MSE at fitted ncomp
%       .cv_mse0        - Intercept-only CV MSE
%       .n_cells        - Cells retained after outlier removal
%       .parms          - Feature names
%       .perm           - Permutation results (if scramble=true):
%           .R2_dist    [1 x n_perm]  - Scrambled R² distribution
%           .cv_mse_dist[1 x n_perm]  - Scrambled CV MSE distribution
%           .p_value    [scalar]       - Fraction of scrambled R² >= real R²

% --- Parse options ---
op.ncomp_results = [];
op.ncomp         = [];
op.cv            = 10;
op.cell_lines    = {'AEN','NEJ'};
op.scramble      = true;
op.n_perm        = 100;
op.use_1se       = false;
op.rm_xsout      = true;  % Pass to pls.m — disable T² outlier removal if false
op = ct_input(varargin, op);

cls = op.cell_lines;
results = struct();

for ci = 1:numel(cls)
    cl = cls{ci};

    % --- Resolve ncomp ---
    nc = resolve_ncomp(op, cl);
    fprintf('\n=== %s: fitting PLS with ncomp=%d, cv=%d ===\n', cl, nc, op.cv);

    % --- Extract data ---
    X     = S.pls.(cl).X{:, S.pls.(cl).parms};
    y     = S.pls.(cl).y;
    parms = S.pls.(cl).parms;
    fprintf('  Cells: %d | Features: %d\n', size(X,1), size(X,2));

    % --- Fit real model ---
    z = pls(X, y, 'ncomp', nc, 'cv', op.cv, 'ploton', false, ...
            'params', parms, 'append', true, 'rm_xsout', op.rm_xsout);

    R2      = sum(z.PCTVAR(2,:));
    cv_mse  = z.MSE(2, end);
    cv_mse0 = z.MSE(2, 1);
    n_out   = size(z.X, 1);

    fprintf('  R² = %.4f | CV MSE = %.4f | n = %d\n', R2, cv_mse, n_out);

    % --- Store real fit ---
    results.(cl).pls_out = z;
    results.(cl).ncomp   = nc;
    results.(cl).R2      = R2;
    results.(cl).cv_mse  = cv_mse;
    results.(cl).cv_mse0 = cv_mse0;
    results.(cl).n_cells = n_out;
    results.(cl).parms   = parms;

    % --- Permutation test ---
    if op.scramble
        fprintf('  Running %d permutations...', op.n_perm);
        R2_dist     = nan(1, op.n_perm);
        cv_mse_dist = nan(1, op.n_perm);

        for pi = 1:op.n_perm
            y_shuf = y(randperm(size(y,1)));  % Shuffle Y only
            z_shuf = pls(X, y_shuf, 'ncomp', nc, 'cv', op.cv, ...
                         'ploton', false, 'params', parms, 'append', true, ...
                         'rm_xsout', op.rm_xsout);
            R2_dist(pi)     = sum(z_shuf.PCTVAR(2,:));
            cv_mse_dist(pi) = z_shuf.MSE(2, end);
        end

        p_value = mean(R2_dist >= R2);
        fprintf(' done. Permutation p = %.4f\n', p_value);

        results.(cl).perm.R2_dist     = R2_dist;
        results.(cl).perm.cv_mse_dist = cv_mse_dist;
        results.(cl).perm.p_value     = p_value;
    end
end

fprintf('\n--- Summary ---\n');
for ci = 1:numel(cls)
    cl = cls{ci};
    r = results.(cl);
    fprintf('%s: ncomp=%d, R²=%.4f, CV MSE=%.4f', cl, r.ncomp, r.R2, r.cv_mse);
    if op.scramble
        fprintf(', perm p=%.4f', r.perm.p_value);
    end
    fprintf('\n');
end

end


%% ========================================================================
function nc = resolve_ncomp(op, cl)
% Resolve ncomp from manual override or ncomp_results struct.
    if ~isempty(op.ncomp)
        nc = op.ncomp;
    elseif ~isempty(op.ncomp_results)
        if op.use_1se
            nc = op.ncomp_results.(cl).ncomp_1se;
        else
            nc = op.ncomp_results.(cl).ncomp_opt;
        end
    else
        error('fit_pls_model:no_ncomp', ...
              'Provide either ''ncomp'' or ''ncomp_results''. No silent defaults.');
    end
end