function Z = eval_surface(xg, yg, fit)
Z = fit.b(1) + fit.b(2).*xg + fit.b(3).*yg;
idx = 4;
for t = 4:numel(fit.terms)
    switch fit.terms{t}
        case 'x2'
            Z = Z + fit.b(t).*xg.^2;
        case 'y2'
            Z = Z + fit.b(t).*yg.^2;
        case 'xy'
            Z = Z + fit.b(t).*(xg.*yg);
    end
end
end