function [minxs, maxxs] = getExtremumsUniform(f, a, b, delta)
    minxs = [];
    maxxs = [];
    
    last = a;
    while(last <= b)
       [extrx, ismin, success] = getExtremumUniform(f, last, b, delta);
       
       if(~success)
           break;
       end
       
       if(ismin == 1)
           minxs = [minxs extrx];
           last = extrx;
       else
           maxxs = [maxxs extrx];
           last = extrx;
       end
    end
end

function [extrx, ismin, success] = getExtremumUniform(f, a, b, delta)
    ismin = 1;
    extrx = b;
    success = 0;
    
    xs = a:delta:b;
    breakCondition = @(i)(f(xs(i + 1)) >= f(xs(i)));
    
    if(f(a) < f(a + delta))
        ismin = 0;
        breakCondition = @(i)(f(xs(i + 1)) <= f(xs(i)));
    end
    
    for i = 1:(length(xs) - 1)
        if(breakCondition(i))
            extrx = xs(i) + delta / 2;
            success = 1;
            break;
        end
    end
end

