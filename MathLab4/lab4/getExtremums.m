function [minxs, maxxs] = getExtremums(f, a, b)
    fwrapmin = @(x)(double(f(x)));
    fwrapmax = @(x)(double(-f(x)));
    
    minxs = []; maxxs = [];
    
    midpointsEq = diff(f, 2) == 0;
    midpoints = solve(midpointsEq);
    
    if(f(midpoints(1) - eps) >= f(midpoints(1)))
        maxxs = [maxxs appendFoundUnimodalExtremum(fwrapmax, a, midpoints(1))];
    else
        minxs = [minxs appendFoundUnimodalExtremum(fwrapmin, a, midpoints(1))];
    end
    
    for i = 1:(length(midpoints) - 1)
        if(f(midpoints(i) + eps) >= f(midpoints(i)))
            maxxs = [maxxs appendFoundUnimodalExtremum(fwrapmax, midpoints(i), midpoints(i + 1))];
        else
            minxs = [minxs appendFoundUnimodalExtremum(fwrapmin, midpoints(i), midpoints(i + 1))];
        end
    end
    
    if(f(midpoints(length(midpoints)) + eps) >= f(midpoints(length(midpoints))))
        maxxs = [maxxs appendFoundUnimodalExtremum(fwrapmax, midpoints(length(midpoints)), b)];
    else
        minxs = [minxs appendFoundUnimodalExtremum(fwrapmin, midpoints(length(midpoints)), b)];
    end
end

function [extrx] = appendFoundUnimodalExtremum(fwrap, a, b)
    extrx = [];
    
    if(b >= a)
        xs = fmincon(fwrap, double(a), [], [], [], [], double(a), double(b));
        if(xs >= a && xs <= b)
            extrx = xs;
        end
    end
end

