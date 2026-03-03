function res = plsr_cv_predvar(T, cellLine, cellLineVar, yVar, xVars, varargin)
% PLSR_CV_PREDVAR  Predicted variance explained (Q^2) by k-fold CV for a PLSR model.
%
% res = plsr_cv_predvar(T, "AEN", "cell_line", "IL8c_log", ["Mean_nfkb_n",...])
%
% NAME-VALUE
%   'K'        : number of folds (default 10)
%   'ncomp'    : number of components (default min(10, numel(xVars)))
%   'seed'     : RNG seed (default 1)
%
% OUTPUT fields:
%   res.pve_pred  : predicted variance explained (Q^2)
%   res.sse_pred  : sum squared prediction error
%   res.sst       : total sum squares
%   res.n         : number of samples used
%   res.ncomp     : n components used
%   res.xVars     : predictors used

p.K = 10;
p.ncomp = [];
p.seed = 1;
p = ct_input(varargin, p);

cellLine = string(cellLine);
cellLineVar = string(cellLineVar);
yVar = string(yVar);
xVars = string(xVars);

Tc = T(T.(cellLineVar) == cellLine, :);

X = Tc{:, xVars};
y = Tc.(yVar);

% Drop NaNs
good = ~any(isnan(X),2) & ~isnan(y);
X = X(good,:); y = y(good);

n = size(X,1);
assert(n > p.K, 'Too few samples (%d) for K=%d CV.', n, p.K);

% Choose components
if isempty(p.ncomp)
    ncomp = min(10, size(X,2));
else
    ncomp = min(p.ncomp, size(X,2));
end

% Create folds
rng(p.seed);
fold = crossvalind('Kfold', n, p.K);

yhat = nan(n,1);

for k = 1:p.K
    te = (fold == k);
    tr = ~te;

    Xtr = X(tr,:); ytr = y(tr);
    Xte = X(te,:);

    % Standardize using training stats
    [Xtrz, xm, xs] = zscore(Xtr);
    ytrm = mean(ytr); ytrs = std(ytr);
    if ytrs == 0, ytrs = 1; end
    ytrz = (ytr - ytrm) ./ ytrs;

    Xtez = (Xte - xm) ./ xs;

    % Fit PLS on training
    [~,~,~,~,beta] = plsregress(Xtrz, ytrz, ncomp);

    % Predict on test (in z-space), then unscale
    ytez = [ones(size(Xtez,1),1), Xtez] * beta;
    yhat(te) = ytez * ytrs + ytrm;
end

sse = nansum((y - yhat).^2);
sst = nansum((y - mean(y)).^2);

pve_pred = 1 - sse/sst;

res = struct('pve_pred', pve_pred, 'sse_pred', sse, 'sst', sst, ...
             'n', n, 'ncomp', ncomp, 'xVars', xVars, 'cellLine', cellLine);
end