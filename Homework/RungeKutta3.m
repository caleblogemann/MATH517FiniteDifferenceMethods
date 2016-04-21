function [u, h, k] = RungeKutta3(a, N, T, icFun)
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
    u(:,1) = icFun(x(0:N+1));

    f = @(u) uxFun(u, a, h);

    for n = 1:Nt
        Y1 = f(u(:,n));
        Y2 = f(u(:,n) + k/2*Y1);
        Y3 = f(u(:,n) + 3*k/4*Y2);
        u(:,n+1) = u(:,n) + k/9*(2*Y1 + 3*Y2 + 4*Y3);
        %plot(x(0:N+1),u(:,n+1));
        %title(num2str(N));
        %pause(k)
    end
end
