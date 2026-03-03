function Tcond = build_condition_means(S, meanCols, zVar)
% BUILD_CONDITION_MEANS  Aggregate per-cell data into condition means.
%
% Tcond = build_condition_means(S, meanCols, zVar)
%   meanCols: cellstr of columns from S.scatter.X (e.g. {'Mean_ekar',...})
%   zVar:     output column from S.Y (e.g. 'IL8c_log')

Tcell = [S.meta, S.scatter.X, S.Y];

Tcond = groupsummary(Tcell, {'cell_line','tx_label'}, 'mean', [meanCols, {zVar}]);

% groupsummary prefixes numeric variables with 'mean_'
Tcond.Properties.VariableNames = strrep(Tcond.Properties.VariableNames, 'mean_', '');
end