% BUILD_META_FROM_P  Construct per-cell metadata table from pulse objects.
%
%   meta = BUILD_META_FROM_P(p)
%
% INPUT
%   p : 1×N cell array of pulse structs (output of build_intermediate_data)
%       Each p{i} must contain:
%           - p{i}.tx   : [nCells × 2] cell array
%                         col 1 = treatment label (char/cellstr)
%                         col 2 = cell line label ('AEN' or 'NEJ')
%           - p{i}.cid  : [nCells × 3] numeric cell identifiers
%           - p{i}.pulse_mat : [nCells × nFeatures] predictor matrix
%
% OUTPUT
%   meta : table [Ntotal × 6]
%       Variables:
%           - exp_id    : experiment index (1–N)
%           - tx_label  : treatment label (string)
%           - cell_line : cell line label (string)
%           - cid1      : numeric cell identifier column 1
%           - cid2      : numeric cell identifier column 2
%           - cid3      : numeric cell identifier column 3
%
% NOTES
%   - One row per cell.
%   - Row ordering matches the concatenated order of p{i}.pulse_mat.
%   - Asserts internal consistency of row counts.
%
% See also: extract_features_from_p, build_il8_from_d
%

function meta = build_meta_from_p(p)
meta = cell(numel(p),1);

for i = 1:numel(p)
    n = size(p{i}.tx,1);

    assert(size(p{i}.cid,1) == n, 'cid/tx mismatch exp %d', i);
    assert(size(p{i}.pulse_mat,1) == n, 'pulse_mat/tx mismatch exp %d', i);

    meta{i} = table( ...
        repmat(i,n,1), ...
        string(p{i}.tx(:,1)), ...
        string(p{i}.tx(:,2)), ...
        p{i}.cid(:,1), p{i}.cid(:,2), p{i}.cid(:,3), ...
        'VariableNames', {'exp_id','tx_label','cell_line','cid1','cid2','cid3'} );
end

meta = vertcat(meta{:});
end

