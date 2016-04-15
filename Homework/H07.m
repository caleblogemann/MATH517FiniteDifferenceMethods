%% Problem 1
    uExact = @(t, x, y) exp(-32*pi^2*t)*cos(4*pi*x).*cos(4*pi*y);
    f = @(x, y) uExact(0, x, y);
    T = 1;
    E = [];
    H = [];
    for N = 10*2.^(1:5) - 1
        [u, Ux, Uy, k] = ADI(N, 1, f, uExact);
        % create exact solution
        t = @(n) n*k;
        uExactMatrix = cell2mat(arrayfun(@(n) uExact(t(n), Ux, Uy)', 0:T/k, 'UniformOutput', false));
        H = [H; Uy(2) - Uy(1)];
        E = [E; norm(u(:,end) - uExactMatrix(:,end), inf)];
    end
    hRatios = H(1:end-1)./H(2:end);
    errorRatios = E(1:end-1)./E(2:end);
    order = log(errorRatios)./log(hRatios);
    table(hRatios, errorRatios, order)
