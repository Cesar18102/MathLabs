 function [xs, ys] = linearInterpolation(xp, yp, delta)
    xs = xp(1):delta:xp(length(xp));
    ys = zeros(size(xs));
    
    ys(1) = yp(1);
    ys(length(ys)) = yp(length(yp));
    
    for i = 2:(length(xs) - 1)
        ys(i) = linearInterpolationStep(xp, yp, xs(i));
    end
end

function [y] = linearInterpolationStep(xp, yp, x)
    i = findIntrpolationIndex(xp, x);
    y = yp(i) + (yp(i + 1) - yp(i)) * (x - xp(i)) / (xp(i + 1) - xp(i));
end

function [i] = findIntrpolationIndex(xp, x)
    for i = 1:(length(xp) - 1)
        if(xp(i) <= x && xp(i + 1) >= x)
            break;
        end
    end
end