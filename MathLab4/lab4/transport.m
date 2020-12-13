function [transfers] = transport(deliveryPrices, storage, demand)
    strgc = height(storage);

    dp = [ deliveryPrices, zeros(strgc, 1) ];
    dmnd = [ demand, (sum(storage) - sum(demand)) ];
    dpv = reshape(dp', 1, []);
    
    argc = width(dpv);
    dmndc = width(dmnd);
    
    f = @(x)(dpv * x);
    
    x0 = ones(argc, 1);
    Aeq = [];
    beq = [];
    
    for i = 1:strgc
        aeq = [ zeros(1, (i - 1) * dmndc), ones(1, dmndc), zeros(1, argc - i * dmndc) ];
        
        Aeq = [ Aeq; aeq ];
        beq = [ beq; storage(i) ];
    end
    
    for i = 1:dmndc
        aeq = [ zeros(strgc, i - 1), ones(strgc, 1), zeros(strgc, dmndc - i) ];
        
        Aeq = [ Aeq; reshape(aeq', 1, []) ];
        beq = [ beq; dmnd(i) ];
    end
    
    temp = fmincon(f, x0, [], [], Aeq, beq, zeros(argc, 1), Inf(argc, 1));
    transfers = reshape(temp, size(dp'))';
end


