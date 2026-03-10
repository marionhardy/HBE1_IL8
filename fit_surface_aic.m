function fit = fit_surface_aic(x, y, z, varargin)
% FIT_SURFACE_AIC  Select polynomial surface (order<=2) by AIC.
%
%   Tests 4 candidate polynomial surfaces (matching legacy pipeline):
%     poly11: x + y           (linear)
%     poly21: x + y + x²      (quadratic in x)
%     poly12: x + y + y²      (quadratic in y)
%     poly22: x + y + x² + y² + xy  (full quadratic)
%
%   fit = fit_surface_aic(x, y, z)
%   fit = fit_surface_aic(x, y, z, 'force', forceSpec)
%
%   forceSpec:
%     "" (default): AIC selection among candidate models
%     "poly11", "poly12", "poly21", "poly22": force specific model
%
% OUTPUTS
%   fit.b        - Coefficients
%   fit.terms    - Term labels
%   fit.SSE      - Sum of squared errors
%   fit.R2       - R²
%   fit.adjR2    - Adjusted R²
%   fit.AIC      - AIC value
%   fit.model    - Model name (poly11, poly12, poly21, poly22)
%   fit.inc_x2   - Includes x² term
%   fit.inc_y2   - Includes y² term
%   fit.inc_xy   - Includes xy term

p.force = "";
p = ct_input(varargin, p);

% 4 candidate models matching legacy pipeline
% Each row: [inc_x2, inc_y2, inc_xy], model name
models = {
    0, 0, 0, 'poly11'   % x + y
    1, 0, 0, 'poly21'   % x + y + x²
    0, 1, 0, 'poly12'   % x + y + y²
    1, 1, 1, 'poly22'   % x + y + x² + y² + xy
};

% Force specific model if requested
if strlength(p.force) > 0
    match = strcmp(models(:,4), p.force);
    if ~any(match)
        error('Unknown force spec: %s. Use poly11, poly12, poly21, or poly22.', p.force);
    end
    models = models(match, :);
end

n = numel(z);
best = [];

for i = 1:size(models, 1)
    inc_x2 = models{i,1};
    inc_y2 = models{i,2};
    inc_xy = models{i,3};
    mname  = models{i,4};

    % Build design matrix
    X_mat = [ones(size(x)), x, y];
    terms = {'1','x','y'};
    if inc_x2; X_mat = [X_mat, x.^2]; terms{end+1} = 'x2'; end
    if inc_y2; X_mat = [X_mat, y.^2]; terms{end+1} = 'y2'; end
    if inc_xy; X_mat = [X_mat, x.*y]; terms{end+1} = 'xy'; end

    b    = X_mat \ z;
    zhat = X_mat * b;
    SSE  = sum((z - zhat).^2);
    SST  = sum((z - mean(z)).^2);
    R2   = 1 - SSE/SST;
    k    = numel(b);
    adjR2 = 1 - (SSE/SST) * (n - 1) / (n - k);  % MATLAB fit() convention
    AIC  = n * log(SSE/n) + 2 * k;

    cur = struct('b', b, 'terms', {terms}, 'SSE', SSE, ...
                 'R2', R2, 'adjR2', adjR2, 'AIC', AIC, ...
                 'model', mname, ...
                 'inc_x2', inc_x2, 'inc_y2', inc_y2, 'inc_xy', inc_xy);

    if isempty(best) || cur.AIC < best.AIC
        best = cur;
    end
end

fit = best;
end