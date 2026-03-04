function [screen] = run_sensor_screen(S, varargin)
% RUN_SENSOR_SCREEN  Combinatorial sensor screen for PLS models.
%
%   Fits PLS for every combination of sensors (1 through all) per cell
%   line. Reports R² and CV MSE for each. Optionally runs permutation
%   test on selected models.
%
% INPUTS
%   S               - Main analysis struct with S.pls.(cl).X, .parms, .y
%
% OPTIONAL (name-value via ct_input)
%   ncomp           - Components to use (default: 10)
%   cv              - CV folds (default: 10)
%   cell_lines      - Cell lines to run (default: {'AEN','NEJ'})
%   perm_top        - Number of top models (by CV MSE) to run permutation
%                     test on (default: 3)
%   n_perm          - Permutations per model (default: 100)
%   plot            - Generate summary plot (default: true)
%
% OUTPUTS
%   screen          - Struct per cell line, each containing a table:
%                     sensors, n_features, R2, cv_mse, n_cells, perm_p

% --- Sensor map (defaults) ---
sensor_default.AEN = {'ampk','erkt','nfkb','nnuc','nnfc'};
sensor_default.NEJ = {'jnkt','ekar','nfkb','nnuc','nnfc'};

% --- Parse options ---
op.ncomp      = 10;
op.cv         = 10;
op.cell_lines = {'AEN','NEJ'};
op.perm_top   = 3;
op.n_perm     = 100;
op.plot       = true;
op.sensors    = [];  % Override sensor list (applies to all cell lines)
op = ct_input(varargin, op);

cls = op.cell_lines;
screen = struct();

for ci = 1:numel(cls)
    cl = cls{ci};
    if ~isempty(op.sensors)
        sensors = op.sensors;
    else
        sensors = sensor_default.(cl);
    end
    nsens   = numel(sensors);

    % --- Generate all combinations (1 through nsens) ---
    combos = {};
    for k = 1:nsens
        c = nchoosek(1:nsens, k);
        for ri = 1:size(c,1)
            combos{end+1} = sensors(c(ri,:)); %#ok<AGROW>
        end
    end
    ncombo = numel(combos);
    fprintf('\n=== %s: %d sensor combinations ===\n', cl, ncombo);

    % --- Extract sensor tags from parms ---
    parms    = S.pls.(cl).parms;
    all_tags = regexp(parms, '_([a-z]+)$', 'tokens');
    all_tags = cellfun(@(x) x{1}{1}, all_tags, 'UniformOutput', false);

    % --- Full X matrix (numeric) ---
    X_full = S.pls.(cl).X{:, parms};
    y      = S.pls.(cl).y;

    % --- Preallocate results ---
    combo_labels = cell(ncombo, 1);
    n_features   = zeros(ncombo, 1);
    R2           = zeros(ncombo, 1);
    cv_mse       = zeros(ncombo, 1);
    n_cells      = zeros(ncombo, 1);

    % --- Sweep combinations ---
    for mi = 1:ncombo
        combo = combos{mi};
        combo_labels{mi} = strjoin(combo, '+');

        % Select columns matching this sensor set
        col_idx    = ismember(all_tags, combo);
        X_sub      = X_full(:, col_idx);
        parms_sub  = parms(col_idx);
        n_features(mi) = sum(col_idx);

        % Cap ncomp at n_features
        nc = min(op.ncomp, n_features(mi));

        % Fit
        z = pls(X_sub, y, 'ncomp', nc, 'cv', op.cv, 'ploton', false, ...
                'params', parms_sub, 'append', true);

        R2(mi)      = sum(z.PCTVAR(2,:));
        cv_mse(mi)  = z.MSE(2, end);
        n_cells(mi) = size(z.X, 1);

        fprintf('  [%2d/%2d] %-30s  nf=%2d  nc=%2d  R²=%.4f  CVE=%.4f  n=%d\n', ...
                mi, ncombo, combo_labels{mi}, n_features(mi), nc, ...
                R2(mi), cv_mse(mi), n_cells(mi));
    end

    % --- Build results table ---
    T = table(combo_labels, n_features, R2, cv_mse, n_cells, ...
              'VariableNames', {'sensors','n_features','R2','cv_mse','n_cells'});
    T = sortrows(T, 'cv_mse', 'ascend');

    % --- Permutation test on top models ---
    perm_p = nan(height(T), 1);
    n_top  = min(op.perm_top, height(T));
    fprintf('  Running permutation test on top %d models...\n', n_top);

    for ti = 1:n_top
        combo    = strsplit(T.sensors{ti}, '+');
        col_idx  = ismember(all_tags, combo);
        X_sub    = X_full(:, col_idx);
        parms_sub = parms(col_idx);
        nc       = min(op.ncomp, sum(col_idx));
        real_R2  = T.R2(ti);

        R2_shuf = nan(1, op.n_perm);
        for pi = 1:op.n_perm
            y_shuf = y(randperm(numel(y)));
            z_shuf = pls(X_sub, y_shuf, 'ncomp', nc, 'cv', op.cv, ...
                         'ploton', false, 'params', parms_sub, 'append', true);
            R2_shuf(pi) = sum(z_shuf.PCTVAR(2,:));
        end
        perm_p(ti) = mean(R2_shuf >= real_R2);
        fprintf('    %-30s  perm p=%.4f\n', T.sensors{ti}, perm_p(ti));
    end
    T.perm_p = perm_p;

    screen.(cl) = T;

    % --- Print summary ---
    fprintf('\n  %s — Top 5 models by CV MSE:\n', cl);
    disp(T(1:min(5,height(T)), :));
end

% --- Plot ---
if op.plot
    plot_sensor_screen(screen, cls);
end

end