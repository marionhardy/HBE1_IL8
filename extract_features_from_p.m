% EXTRACT_FEATURES_FROM_P  Extract selected predictor columns from pulse data.
%
%   F = EXTRACT_FEATURES_FROM_P(p, feature_names)
%
% INPUT
%   p : 1×N cell array of pulse structs
%       Required fields:
%           - p{i}.pulse_mat : [nCells × nFeatures]
%           - p{i}.parm      : 1×nFeatures cell array of feature names
%
%   feature_names : cellstr
%       Exact feature names to extract (e.g., {'Mean_ekar','Mean_ampk'})
%
% OUTPUT
%   F : table [Ntotal × numel(feature_names)]
%       Each column corresponds to a requested feature.
%       Rows are concatenated across experiments in p.
%
% NOTES
%   - Missing features in an experiment are skipped with warning.
%   - Row ordering matches build_meta_from_p(p).
%   - Variable names are sanitized using matlab.lang.makeValidName.
%
% Example
%   feats = {'Mean_ekar','Mean_ampk'};
%   F = extract_features_from_p(p, feats);
%
% See also: build_meta_from_p
%

function F = extract_features_from_p(p, feature_names)
Fcell = cell(numel(p),1);

for i = 1:numel(p)
    X = p{i}.pulse_mat;
    parms = string(p{i}.parm);

    [tf, loc] = ismember(string(feature_names), parms);
    missing = feature_names(~tf);
    if ~isempty(missing)
        warning('Exp %d missing %d features (e.g. %s)', i, numel(missing), missing{1});
    end

    loc = loc(tf);
    fn  = matlab.lang.makeValidName(string(feature_names(tf)));

    Fi = array2table(X(:,loc), 'VariableNames', cellstr(fn));
    Fi.exp_id = repmat(i, size(X,1), 1);

    Fcell{i} = Fi;
end

F = vertcat(Fcell{:});
end
