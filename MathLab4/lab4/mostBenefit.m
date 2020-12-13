function [benefits] = mostBenefit(timeSpend, timeLeft, profit)
    argc = width(profit);
    f = @(x)(-profit * x);
    x0 = zeros(argc, 1);
    benefits = fmincon(f, x0, timeSpend, timeLeft, [], [], zeros(argc, 1), Inf(argc, 1));
end

