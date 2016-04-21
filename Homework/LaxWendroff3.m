function [u, h, k] = LaxWendroff3(a, N, T, f)
    % discretize space
    h = 1/(N+1);
    x = @(i) i*h;

    % discretize time
    % find k such that k = O(h), but not exact k = h
    Nt = ceil(T/(.9*h));
    k = T/Nt;
    t = @(n) n*k;

    % initiate matrix to store u at times 1:Nt
    u = zeros(N+2, Nt+1);
    % add initial conditions
    u(:,1) = f(x(0:N+1));

    for n=1:Nt
        % create U^n_{j-2}, U^n_{j-1}, and U^n_{j+1} considering periodic
        % boundary conditions
        Uj2m = [u(end-1,n); u(end, n); u(1:end-2,n)];
        Uj1m = [u(end, n); u(1:end-1,n)];
        Uj1p = [u(2:end,n); u(1,n)];
        u(:,n+1) = u(:,n) - (a*k)/(6*h)*(Uj2m - 6*Uj1m + 3*u(:,n) + 2*Uj1p)...
            + (a*k)^2/(2*h^2)*(Uj1m - 2*u(:,n) + Uj1p)...
            + (a*k)^3/(6*h^3)*(Uj2m - 3*Uj1m + 3*u(:,n) - Uj1p);
    end
end
