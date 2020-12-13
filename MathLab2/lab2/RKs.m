function [res] = RKs(fs, a, x0, y0, h)
    xv = (x0:h:a)';
    
    yv = zeros(length(xv), numel(y0));
    yv(1,:) = y0(:,1); 
    
    for i=1:(length(xv)-1)
        ys1 = num2cell(yv(i,:));
        k1 = fs(ys1{:})'; 
        
        ys2 = num2cell(yv(i,:) + k1 * 0.5 * h);
        k2 = fs(ys2{:})'; 
        
        ys3 = num2cell(yv(i,:) + k2 * 0.5 * h);
        k3 = fs(ys3{:})'; 
        
        ys4 = num2cell(yv(i,:) + k3 * h);
        k4 = fs(ys4{:})'; 
        
        yv(i + 1,:) = yv(i,:) + (k1 + 2 * k2 + 2 * k3 + k4) * (h / 6); 
    end
    
    res = [xv, yv];
end 