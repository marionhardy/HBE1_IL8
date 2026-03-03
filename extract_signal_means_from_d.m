function Xtbl = extract_signal_means_from_d(d, signals, varargin)
% EXTRACT_SIGNAL_MEANS_FROM_D  Compute per-cell signal summaries from decompressed d.
%
%   Xtbl = EXTRACT_SIGNAL_MEANS_FROM_D(d, signals, 'tidx', tidx)
%
% INPUT
%   d       : 1×N cell array; each d{i} is a compressed dataset consumable by ct_compform
%   signals : cellstr of channel names in pinfo.cname (e.g., {'ekar','erktr','jnktr','ampkar','nfkb_n_fc'})
%
% NAME-VALUE
%   tidx    : vector of time indices to average over (default = all timepoints)
%   fun     : summary function for time dimension: 'mean' (default) or 'median' or 'max'
%
% OUTPUT
%   Xtbl : table [Ntotal × (numel(signals)+1)]
%       Columns:
%         - exp_id
%         - Mean_<signal> (or Median_/Max_ depending on fun)
%
% NOTES
%   - Uses ct_compform(d{i}), then concatenates blocks with vertcat(dt_ar{:}).
%   - Assumes dt blocks are [cells × time × channels].
%
p.tidx = [];
p.fun  = 'mean';
p = ct_input(varargin, p);

Xcell = cell(numel(d),1);

for i = 1:numel(d)
    [dt_ar, ~, pinfo] = ct_compform(d{i});
    D = vertcat(dt_ar{:});                % [nCells × nTime × nCh]
    [nCells, nTime, ~] = size(D);

    if isempty(p.tidx)
        tidx = 1:nTime;
    else
        tidx = p.tidx(:)';
        tidx = tidx(tidx>=1 & tidx<=nTime);
        assert(~isempty(tidx), 'tidx invalid for exp %d', i);
    end

    cname = string(pinfo.cname);
    Ti = table();

    for s = 1:numel(signals)
        sig = string(signals{s});
        ch = find(cname == sig, 1);
        assert(~isempty(ch), 'Signal "%s" not found in exp %d', sig, i);

        X = D(:, tidx, ch);  % [nCells × numel(tidx)]

        switch lower(p.fun)
            case 'mean'
                v = mean(X, 2, 'omitnan');
                vname = "Mean_" + sig;
            case 'median'
                v = median(X, 2, 'omitnan');
                vname = "Median_" + sig;
            case 'max'
                v = max(X, [], 2);
                vname = "Max_" + sig;
            otherwise
                error('Unknown fun: %s', p.fun);
        end

        Ti.(matlab.lang.makeValidName(vname)) = v;
    end

    Xcell{i} = Ti;
end

Xtbl = vertcat(Xcell{:});
end