%% Problem 1
    uExact = @(t, x, y) exp(-32*pi^2*t)*cos(4*pi*x).*cos(4*pi*y);
    f = @(x, y) uExact(0, x, y);
    T = 1;
    E = [];
    H = [];
    for N = 10*2.^(6:10) - 1
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

%% Problem 4 for LaxWendroff3 created in Problem 2
    a = 1;
    T = 1;
    uExact = @(t, x) 2*exp(-200*(mod((x-a*t),1) - 1/2).^2);
    f = @(x) uExact(0,x);
    E = [];
    H = [];
    for N = 10*2.^(1:8) - 1
        [u, h, k] = LaxWendroff3(a, N, T, f);
        % create exact solution
        t = @(n) n*k;
        x = @(i) i*h;
        uExactMatrix = cell2mat(arrayfun(@(n) uExact(t(n), x(0:N+1))', 0:T/k, 'UniformOutput', false));
        H = [H; h];
        E = [E; norm(u(:,end) - uExactMatrix(:,end), inf)];
    end
    hRatios = H(1:end-1)./H(2:end);
    errorRatios = E(1:end-1)./E(2:end);
    order = log(errorRatios)./log(hRatios);
    table(hRatios, errorRatios, order)

    uExact = @(t, x) 2*exp(-200*(mod((x-a*t),1) - 1/2).^2).*cos(40*pi*mod((x-a*t),1));
    f = @(x) uExact(0,x);
    E = [];
    H = [];
    for N = 10*2.^(1:8) - 1
        [u, h, k] = LaxWendroff3(a, N, T, f);
        % create exact solution
        t = @(n) n*k;
        x = @(i) i*h;
        uExactMatrix = cell2mat(arrayfun(@(n) uExact(t(n), x(0:N+1))', 0:T/k, 'UniformOutput', false));
        H = [H; h];
        E = [E; norm(u(:,end) - uExactMatrix(:,end), inf)];
    end
    hRatios = H(1:end-1)./H(2:end);
    errorRatios = E(1:end-1)./E(2:end);
    order = log(errorRatios)./log(hRatios);
    table(hRatios, errorRatios, order)

%% Problem 4 for RungeKutta3 created in Problem 3
    a = 1;
    T = 1;
    uExact = @(t, x) 2*exp(-200*(mod((x-a*t),1) - 1/2).^2);
    f = @(x) uExact(0,x);
    E = [];
    H = [];
    for N = 10*2.^(1:8) - 1
        [u, h, k] = RungeKutta3(a, N, T, f);
        % create exact solution
        t = @(n) n*k;
        x = @(i) i*h;
        uExactMatrix = cell2mat(arrayfun(@(n) uExact(t(n), x(0:N+1))', 0:T/k, 'UniformOutput', false));
        H = [H; h];
        E = [E; norm(u(:,end) - uExactMatrix(:,end), inf)];
    end
    hRatios = H(1:end-1)./H(2:end);
    errorRatios = E(1:end-1)./E(2:end);
    order = log(errorRatios)./log(hRatios);
    table(hRatios, errorRatios, order)

    uExact = @(t, x) 2*exp(-200*(mod((x-a*t),1) - 1/2).^2).*cos(40*pi*mod((x-a*t),1));
    f = @(x) uExact(0,x);
    E = [];
    H = [];
    for N = 10*2.^(4:10) - 1
        [u, h, k] = RungeKutta3(a, N, T, f);
        % create exact solution
        t = @(n) n*k;
        x = @(i) i*h;
        uExactMatrix = cell2mat(arrayfun(@(n) uExact(t(n), x(0:N+1))', 0:T/k, 'UniformOutput', false));
        H = [H; h];
        E = [E; norm(u(:,end) - uExactMatrix(:,end), inf)];
    end
    hRatios = H(1:end-1)./H(2:end);
    errorRatios = E(1:end-1)./E(2:end);
    order = log(errorRatios)./log(hRatios);
    table(hRatios, errorRatios, order)
