function [y] = RK(f, a, x0, y0, h)
    x = x0:h:a;
    
    y = zeros(size(x));
    y(1) = y0;

    for i=1:(length(x)-1)
        k1 = f(x(i), y(i));
        k2 = f(x(i) + 0.5 * h, y(i) + k1 * 0.5 * h);
        k3 = f(x(i) + 0.5 * h, y(i) + k2 * 0.5 * h);
        k4 = f(x(i) + h, y(i) + k3 * h);
        y(i+1) = y(i) + (k1 + 2 * k2 + 2 * k3 + k4) * (h / 6);
    end
end

