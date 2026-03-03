function R = check_feature_redundancy(S, varargin)
% CHECK_FEATURE_REDUNDANCY  Correlation-based redundancy diagnostic for PLS X.
%
% R = check_feature_redundancy(S, 'thresh', 0.95, 'plot', true)
% R = check_feature_redundancy(S, 'thresh', 0.95, 'cell_lines', {'AEN'})
%
% OUTPUT R.(AEN|NEJ):
%   .C           [nFeatures x nFeatures] correlation matrix
%   .parms       feature names (no zero-variance removal — already pruned in S)
%   .flag_pairs  table of redundant pairs sorted by |r|

p.thresh     = 0.95;
p.plot       = true;
p.cell_lines = {'AEN','NEJ'};
p = ct_input(varargin, p);

sensor_tags = {'ekar','ampk','erkt','nfkb','jnkt','nnuc','nnfc'};

for ci = 1:numel(p.cell_lines)
    cl    = p.cell_lines{ci};
    X     = table2array(S.pls.(cl).X(:, 1:end-1));  % drop exp_id
    parms = string(S.pls.(cl).parms);

    % Correlation matrix
    C = corr(X, 'rows', 'pairwise');11
    R.(cl).C     = C;
    R.(cl).parms = parms;

    % Flag redundant pairs — upper triangle only
    C_upper = abs(C);
    C_upper(tril(true(size(C_upper)))) = 0;
    [ri, ci_] = find(C_upper >= p.thresh);

    if isempty(ri)
        R.(cl).flag_pairs = table();
        fprintf('[%s] No redundant pairs at |r| >= %.2f\n', cl, p.thresh);
    else
        R.(cl).flag_pairs = sortrows( ...
            table(parms(ri)', parms(ci_)', ...
                  C_upper(sub2ind(size(C_upper), ri, ci_)), ...
                  'VariableNames', {'Feature_A','Feature_B','abs_r'}), ...
            'abs_r', 'descend');
        fprintf('[%s] %d redundant pairs at |r| >= %.2f\n', ...
            cl, height(R.(cl).flag_pairs), p.thresh);
    end
end

%% Plots
if ~p.plot; return; end

figure('Name','Feature Redundancy','Position',[100 100 700*numel(p.cell_lines) 600]);
for ci = 1:numel(p.cell_lines)
    cl    = p.cell_lines{ci};
    parms = R.(cl).parms;

    [sorted_idx, tick_pos, tick_labels, boundaries] = sensor_order(parms, sensor_tags);
    C_ord = abs(R.(cl).C(sorted_idx, sorted_idx));

    subplot(1, numel(p.cell_lines), ci);
    imagesc(C_ord); colormap(gca, hot); clim([0 1]); colorbar; axis square;
    hold on;
    arrayfun(@(b) xline(b+0.5,'c-','LineWidth',1.2), boundaries(2:end-1));
    arrayfun(@(b) yline(b+0.5,'c-','LineWidth',1.2), boundaries(2:end-1));
    set(gca, 'XTick', tick_pos, 'XTickLabel', tick_labels, 'FontSize', 11, ...
             'YTick', tick_pos, 'YTickLabel', tick_labels);
    xtickangle(45);
    title(sprintf('%s | redundant pairs: %d', cl, height(R.(cl).flag_pairs)));
end
sgtitle('Feature correlation |r| — ordered by sensor block');

% NFkB zoom
figure('Name','NFkB Redundancy','Position',[100 100 700*numel(p.cell_lines) 500]);
for ci = 1:numel(p.cell_lines)
    cl    = p.cell_lines{ci};
    parms = R.(cl).parms;
    idx   = contains(parms, {'nfkb','nnuc','nnfc'});
    C_n   = abs(R.(cl).C(idx, idx));
    lbl   = replace(parms(idx), '_', ' ');

    subplot(1, numel(p.cell_lines), ci);
    imagesc(C_n); colormap(gca, hot); clim([0 1]); colorbar; axis square;
    set(gca, 'XTick', 1:sum(idx), 'XTickLabel', lbl, 'FontSize', 7, ...
             'YTick', 1:sum(idx), 'YTickLabel', lbl);
    xtickangle(45);
    title(sprintf('%s — NFkB block |r|', cl));
end
sgtitle('NFkB representation redundancy: nfkb vs nnuc vs nnfc');
end

%% Helper
function [sorted_idx, tick_pos, tick_labels, boundaries] = sensor_order(parms, tags)
parms_str  = string(parms);
sorted_idx = [];
boundaries = [0];
tick_labels = {};

for s = 1:numel(tags)
    switch tags{s}
        case 'nfkb'
            idx = find(contains(parms_str,'nfkb') & ...
                      ~contains(parms_str,'nnuc') & ...
                      ~contains(parms_str,'nnfc'));
        otherwise
            idx = find(contains(parms_str, tags{s}));
    end
    if isempty(idx); continue; end  % skip absent sensors silently
    sorted_idx   = [sorted_idx;  idx(:)];
    boundaries(end+1) = boundaries(end) + numel(idx);
    tick_labels{end+1} = tags{s};
end

assert(numel(unique(sorted_idx)) == numel(sorted_idx), ...
    'Duplicate indices — check tag overlap');
assert(numel(sorted_idx) == numel(parms), ...
    'sensor_order: %d/%d features assigned', numel(sorted_idx), numel(parms));

tick_pos   = round((boundaries(1:end-1) + boundaries(2:end)) / 2);
boundaries = boundaries;
end