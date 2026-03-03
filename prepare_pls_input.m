function [X_mat, y, parms] = prepare_pls_input(S, R, cell_line, sensors, varargin)
% PREPARE_PLS_INPUT  Build clean [X, y, parms] for PLS, per cell line.
%
% [X_mat, y, parms] = prepare_pls_input(S, R, cell_line, sensors)
% [X_mat, y, parms] = prepare_pls_input(S, R, cell_line, sensors, 'y_var', 'IL8c_raw')
%
% INPUTS
%   S         : main data struct (from loading script)
%   R         : redundancy struct (from check_feature_redundancy)
%   cell_line : 'AEN' or 'NEJ'
%   sensors   : cellstr of sensor tags to include e.g. {'ampk','erkt','nfkb','nnuc','nnfc'}
%
% NAME-VALUE
%   y_var     : Y column from S.Y to use (default: 'IL8c_raw')
%   exclude   : additional feature substrings to exclude (default: {'_tp','mult'})
%               mirrors previous student's pls_col_idx exclude logic
%
% OUTPUTS
%   X_mat  : [nCells x nFeatures] double, z-scored
%   y      : [nCells x 1] double
%   parms  : [1 x nFeatures] cellstr of feature names

p.y_var   = 'IL8c_raw';
p.exclude = {'_tp', 'mult'};
p = ct_input(varargin, p);

%% 1. Cell line mask
assert(ismember(cell_line, {'AEN','NEJ'}), 'cell_line must be AEN or NEJ');
idx = S.idx.(cell_line);

%% 2. Start from valid (non-zero-variance) features for this cell line
valid_parms = string(R.(cell_line).parms_valid);  % already zero-var pruned

%% 3. Keep only requested sensors
% Use same exact-match logic as check_feature_redundancy to avoid tag overlap
keep = false(size(valid_parms));
for s = 1:numel(sensors)
    tag = sensors{s};
    switch tag
        case 'nfkb'
            keep = keep | (contains(valid_parms,'nfkb') & ...
                          ~contains(valid_parms,'nnuc') & ...
                          ~contains(valid_parms,'nnfc'));
        otherwise
            keep = keep | contains(valid_parms, tag);
    end
end

%% 4. Apply exclusions (e.g. time-windowed '_tp', interaction 'mult')
for e = 1:numel(p.exclude)
    keep = keep & ~contains(valid_parms, p.exclude{e});
end

%% 5. Drop redundant features — keep Feature_A, drop Feature_B
% (Feature_A is retained by convention since flag_pairs is sorted by |r|
%  and Feature_A/B ordering is consistent with parms ordering)
pairs = R.(cell_line).flag_pairs;
if ~isempty(pairs)
    drop_redundant = string(pairs.Feature_B);  % drop the second of each redundant pair
    keep = keep & ~ismember(valid_parms, drop_redundant);
end

assert(any(keep), 'No features remain after filtering — check sensor tags');

selected_parms = valid_parms(keep);
fprintf('[%s] Features selected: %d / %d valid\n', ...
    cell_line, sum(keep), numel(valid_parms));

%% 6. Extract X from full feature table
all_parms  = string(S.pls.parms);
X_full     = table2array(S.pls.X(:, 1:end-1));  % drop exp_id

[~, col_idx] = ismember(selected_parms, all_parms);
assert(all(col_idx > 0), 'Some selected features not found in S.pls.X');

X_mat = X_full(idx, col_idx);
parms = cellstr(selected_parms);

%% 7. Z-score X (pls.m also z-scores internally, but explicit here for transparency)
% This ensures X is inspectable in the same scale that enters plsregress
[X_mat, ~, ~] = zscore(X_mat);

%% 8. Extract y
assert(ismember(p.y_var, S.Y.Properties.VariableNames), ...
    'y_var "%s" not found in S.Y', p.y_var);
y = S.Y.(p.y_var)(idx);

%% 9. Final checks
assert(size(X_mat,1) == numel(y),      'X/y row mismatch');
assert(size(X_mat,2) == numel(parms),  'X/parms column mismatch');

fprintf('[%s] Final X: [%d cells x %d features] | y: %s\n', ...
    cell_line, size(X_mat,1), size(X_mat,2), p.y_var);
end