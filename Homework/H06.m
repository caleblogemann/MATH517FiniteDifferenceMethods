%% Problem 2 (a)
    H = [];
    E = [];
    for m = [11, 23, 47, 95, 191, 383]
        [h, k, err] = heat_CN(m, 4);
        H = [H; h];
        E = [E; err];
    end
    hRatios = H(1:end-1)./H(2:end);
    errorRatios = E(1:end-1)./E(2:end);
    order = log(errorRatios)./log(hRatios);
    table(hRatios, errorRatios, order)

%% Problem 2 (b)
    H = [];
    E = [];
    for m = [11, 23, 47, 95, 191, 383]
        [h, k, err] = heat_trbdf2(m, 4);
        H = [H; h];
        E = [E; err];
    end
    hRatios = H(1:end-1)./H(2:end);
    errorRatios = E(1:end-1)./E(2:end);
    order = log(errorRatios)./log(hRatios);
    table(hRatios, errorRatios, order)

%% Problem 3 (a)
    H = [];
    E = [];
    for m = 38*2.^(0:4)
        [h, k, err] = heat_CN3(m, .1);
        H = [H; h];
        E = [E; err];
    end
    hRatios = H(1:end-1)./H(2:end);
    errorRatios = E(1:end-1)./E(2:end);
    order = log(errorRatios)./log(hRatios);
    table(hRatios, errorRatios, order)

%% Problem 3 (b)
    H = [];
    E = [];
    for m = 38*2.^(0:4)
        [h, k, err] = heat_trbdf23(m, 4);
        H = [H; h];
        E = [E; err];
    end
    hRatios = H(1:end-1)./H(2:end);
    errorRatios = E(1:end-1)./E(2:end);
    order = log(errorRatios)./log(hRatios);
    table(hRatios, errorRatios, order)
