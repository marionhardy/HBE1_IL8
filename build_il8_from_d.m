% BUILD_IL8_FROM_D  Extract IL8 (or other channel) response aligned to pulse rows.
%
%   Y = BUILD_IL8_FROM_D(d, p, ch)
%
% INPUT
%   d  : 1×N cell array of experiment data structs
%   p  : 1×N cell array of pulse structs
%   ch : string or char
%        Channel name to extract (e.g., "IL8c")
%
% OUTPUT
%   Y : table [Ntotal × 2]
%       Variables:
%           - <ch>_raw : raw extracted values
%           - <ch>_log : log10-shifted version
%
% PROCESS
%   - Uses make_pls_output_mhy to extract per-experiment channel values.
%   - Asserts that extracted row counts match p{i}.pulse_mat rows.
%   - Concatenates across experiments.
%   - Applies log10 shift: log10(value + abs(floor(min(value))))
%
% ASSUMPTION
%   Row ordering in make_pls_output_mhy matches p{i}.pulse_mat ordering.
%   (Strict length assertions enforced.)
%
% See also: make_pls_output_mhy
%


function Y = build_il8_from_d(d, p, ch)
y_cell = make_pls_output_mhy(ch, d);
assert(iscell(y_cell) && numel(y_cell)==numel(p), 'Unexpected make_pls_output_mhy output');

y_all = [];
for i = 1:numel(p)
    yi = y_cell{i};
    n  = size(p{i}.pulse_mat,1);
    assert(numel(yi) == n, 'IL8 length mismatch exp %d: y=%d p=%d', i, numel(yi), n);
    y_all = [y_all; yi(:)];
end

var_raw = char(string(ch) + "_raw");
var_log = char(string(ch) + "_log");

Y = table(y_all, 'VariableNames', {var_raw});

v = Y.(var_raw);
shift = abs(floor(min(v)));
Y.(var_log) = log10(v + shift);
end
