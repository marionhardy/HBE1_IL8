function plot_pve_subset_bars(PVE, inc, cellLines, rowLabels)
% PLOT_PVE_SUBSET_BARS  Barplot of predicted variance explained for subset models
% with a "+" inclusion matrix underneath.
%
% PVE: struct with fields per cell line (e.g., PVE.AEN, PVE.NEJ), each 7x1
% inc: (nRows x 7) logical/numeric inclusion matrix (rows=metrics, cols=models)
% cellLines: string array, e.g. ["AEN","NEJ"]
% rowLabels: string array of length nRows (y-axis labels for "+" matrix rows)

cellLines = string(cellLines);
rowLabels = string(rowLabels);

nModels = 7;
nRows = size(inc,1);

% Layout similar to original (two panels)
figure;
tiledlayout(1, numel(cellLines), 'TileSpacing','compact', 'Padding','compact');

for j = 1:numel(cellLines)
    cl = cellLines(j);
    nexttile;

    y = PVE.(matlab.lang.makeValidName(cl));
    bar(y); hold on;

    % Match original aesthetic: no x tick labels (we draw "+" matrix instead)
    xlim([0.5, nModels+0.5]);
    xticks(1:nModels);
    xticklabels(repmat("", nModels, 1));
    box on;
    title(cl + " Cell Line", 'Interpreter','none');

    ylabel('Predicted variance explained');

    % Choose y-limits with room below for "+" matrix
    ymax = max([0.25; y(:)]) * 1.05;     % original looked capped around ~0.25
    ymin = -0.10;                        % space for plus matrix
    ylim([ymin, ymax]);

    % Draw "+" matrix under bars
    % Define vertical positions for the rows under 0
    yRow = linspace(-0.02, -0.085, nRows);

    for r = 1:nRows
        for i = 1:nModels
            if inc(r,i)
                text(i, yRow(r), '+', 'HorizontalAlignment','center', ...
                    'VerticalAlignment','middle', 'FontSize', 12);
            end
        end
        % Row label at left
        text(0.35, yRow(r), rowLabels(r), 'HorizontalAlignment','right', ...
            'VerticalAlignment','middle');
    end

    % Optional: show numeric values above bars (comment out to match original more closely)
    % for i = 1:nModels
    %     text(i, y(i)+0.01, sprintf('%.3f', y(i)), 'HorizontalAlignment','center');
    % end

    hold off;
end
end