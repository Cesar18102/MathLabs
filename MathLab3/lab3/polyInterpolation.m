function [xs, ys] = polyInterpolation(x, y, delta)
    w = getW(x);
    as = linsolve(w, y);
    
    if(delta == 0)
        xs = x;
        ys = polyval(flip(as), x);
    else
        xs = x(1):delta:x(length(x));
        ys = polyval(flip(as), xs);
    end
end

function [w] = getW(x)
    len = length(x);
    w = zeros(len, len + 1);
    
    for i = 1:len
        for j = 1:(len + 1)
            w(i, j) = x(i) ^ (j - 1);
        end
    end
end