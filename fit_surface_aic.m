function fit = fit_surface_aic(x, y, z, varargin)
% FIT_SURFACE_AIC  Select polynomial surface (order<=2) by AIC.
%
% fit = fit_surface_aic(x,y,z,'force',forceSpec)
% forceSpec:
%   "" (default): AIC selection among candidate models
%   "NEJ_ERK_NFKB": force {x, y, y^2} (NFkB squared)
%   "AEN_ERK_NFKB": force {x, y, x^2} (ERK squared)

p.force = "";
p = ct_input(varargin, p);

% Candidate term sets (always include intercept, x, y)
% Each row: [x2 y2 xy]
cand = [
    0 0 0  % linear
    1 0 0  % +x^2
    0 1 0  % +y^2
    0 0 1  % +xy
    1 1 0  % +x^2 +y^2
    1 0 1  % +x^2 +xy
    0 1 1  % +y^2 +xy
    1 1 1  % full quadratic
];

% Forced special-case models
if strlength(p.force)>0
    switch string(p.force)
        case "NEJ_ERK_NFKB"
            cand = [0 1 0]; % x,y,y^2
        case "AEN_ERK_NFKB"
            cand = [1 0 0]; % x,y,x^2
        otherwise
            error('Unknown force spec: %s', p.force);
    end
end

best = [];
for i = 1:size(cand,1)
    inc_x2 = cand(i,1); inc_y2 = cand(i,2); inc_xy = cand(i,3);

    [b, zhat, terms] = fit_once(x,y,z,inc_x2,inc_y2,inc_xy);

    SSE = sum((z - zhat).^2);
    SST = sum((z - mean(z)).^2);
    R2  = 1 - SSE/SST;

    n = numel(z);
    k = numel(b);
    AIC = n*log(SSE/n) + 2*k;

    cur = struct('b',b,'terms',{terms},'SSE',SSE,'R2',R2,'AIC',AIC, ...
                 'inc_x2',inc_x2,'inc_y2',inc_y2,'inc_xy',inc_xy);

    if isempty(best) || cur.AIC < best.AIC
        best = cur;
    end
end

fit = best;
end

function [b, zhat, terms] = fit_once(x,y,z,inc_x2,inc_y2,inc_xy)
X = [ones(size(x)), x, y];
terms = {'1','x','y'};

if inc_x2
    X = [X, x.^2]; terms{end+1} = 'x2';
end
if inc_y2
    X = [X, y.^2]; terms{end+1} = 'y2';
end
if inc_xy
    X = [X, x.*y]; terms{end+1} = 'xy';
end

b = X \ z;
zhat = X*b;
end