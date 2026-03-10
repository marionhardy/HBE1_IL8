function [results] = fit_reduced_pls_legacy(S, varargin)
% FIT_REDUCED_PLS_LEGACY  Reduced PLS using previous student's approach.
%
%   Replicates the legacy pipeline: fit full model, select VIP>1 features,
%   refit reduced model on the same data. This is the "double-dipping"
%   approach — feature selection and evaluation use the same data.
%   Provided for direct comparison against the bootstrap-stability method.
%
% INPUTS
%   S               - Main analysis struct with S.pls.(cl).X, .parms, .y
%
% OPTIONAL (name-value via ct_input)
%   ncomp           - Components for full model (default: 10)
%   cv              - CV folds (default: 10)
%   cell_lines      - Cell lines to run (default: {'AEN','NEJ'})
%   rm_xsout        - Pass to pls.m (default: false)
%   vip_thresh      - VIP threshold for feature selection (default: 1)
%
% OUTPUTS
%   results         - Struct per cell line:
%       .selected_parms - Features with VIP > vip_thresh
%       .vip_scores     [1 x nFeatures] - VIP scores from full model
%       .full_R2        - Full model R²
%       .full_cv_mse    - Full model CV MSE
%       .reduced_R2     - Reduced model R²
%       .reduced_cv_mse - Reduced model CV MSE
%       .full_pls       - Full model pls.m output
%       .reduced_pls    - Reduced model pls.m output
%       .n_cells_full   - Cells retained in full model
%       .n_cells_reduced - Cells retained in reduced model
%       .parms          - All feature names

% --- Parse options ---
op.ncomp      = 10;
op.cv         = 10;
op.cell_lines = {'AEN','NEJ'};
op.rm_xsout   = false;
op.vip_thresh = 1;
op = ct_input(varargin, op);

cls = op.cell_lines;
results = struct();

for ci = 1:numel(cls)
    cl = cls{ci};
    fprintf('\n=== %s: Legacy reduced PLS (VIP > %.1f from full model) ===\n', ...
        cl, op.vip_thresh);

    % --- Extract data ---
    parms = S.pls.(cl).parms;
    X     = S.pls.(cl).X{:, parms};
    y     = S.pls.(cl).y;
    [N, nFeat] = size(X);
    fprintf('  Cells: %d | Features: %d\n', N, nFeat);

    % --- Fit full model ---
    fprintf('  Fitting full model (ncomp=%d, CV=%d)...\n', op.ncomp, op.cv);
    z_full = pls(X, y, 'ncomp', op.ncomp, 'cv', op.cv, 'ploton', false, ...
                 'params', parms, 'append', true, 'rm_xsout', op.rm_xsout);
    full_R2     = sum(z_full.PCTVAR(2,:));
    full_cv_mse = z_full.MSE(2, end);

    % --- VIP selection (single fit, no resampling) ---
    [~, ~, vip_all] = pls_VIP_score(z_full, parms);
    selected_mask  = vip_all' > op.vip_thresh;
    selected_parms = parms(selected_mask);

    fprintf('  Selected: %d / %d features with VIP > %.1f\n', ...
        sum(selected_mask), nFeat, op.vip_thresh);

    % --- Refit reduced model on SAME data (legacy approach) ---
    nc_reduced = min(op.ncomp, sum(selected_mask));
    fprintf('  Refitting reduced model (%d features, ncomp=%d, CV=%d)...\n', ...
        sum(selected_mask), nc_reduced, op.cv);
    z_reduced = pls(X(:, selected_mask), y, 'ncomp', nc_reduced, 'cv', op.cv, ...
                    'ploton', false, 'params', selected_parms, 'append', true, ...
                    'rm_xsout', op.rm_xsout);
    reduced_R2     = sum(z_reduced.PCTVAR(2,:));
    reduced_cv_mse = z_reduced.MSE(2, end);

    % --- Report ---
    fprintf('\n  %s Results (LEGACY — same-data selection):\n', cl);
    fprintf('    Full model:    R²=%.4f  CV MSE=%.4f  n=%d  (%d features)\n', ...
        full_R2, full_cv_mse, size(z_full.X,1), nFeat);
    fprintf('    Reduced model: R²=%.4f  CV MSE=%.4f  n=%d  (%d features)\n', ...
        reduced_R2, reduced_cv_mse, size(z_reduced.X,1), sum(selected_mask));
    fprintf('    NOTE: Feature selection used same data as evaluation (double-dipping).\n');

    % --- Store ---
    results.(cl).selected_parms  = selected_parms;
    results.(cl).selected_mask   = selected_mask;
    results.(cl).vip_scores      = vip_all';
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

end