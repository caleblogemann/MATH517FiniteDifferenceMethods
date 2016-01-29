%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Problem 6
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% solving BVP u'' + u = f(x) for f(x) = -e^x
% with BCs u'(0) - u(0) = 0 and u'(L) + u(L) = 0, where L = 10

L = 10;
f = @(x) -exp(x);

% exact solution
u = @(x) (-1/2)*exp(x) + (exp(L)/(2*cos(L)))*(sin(x) + cos(x));

% create vector for storing errors at respective h values
E = [];
% N is number of interior point
for N=[10, 20, 40, 80, 160, 320, 640, 1280, 2560, 5120]
    % find spacing
    h = L/(N+1);

    % create x grid and function values on grid
    x = 0:h:L;
    fx = f(x(2:N+1))';

    % create matrix A that solves the finite difference problem
    e = ones(N, 1);
    mid = (-2 + h^2)*e;
    mid(1) = mid(1) + 4/(3 + 2*h);
    mid(N) = mid(N) + 4/(3 + 2*h);
    low = e;
    low(end - 1) = 1 - 1/(3 + 2*h);
    up = e;
    up(2) = 1 - 1/(3 + 2*h);
    A = 1/h^2 * spdiags([low, mid, up], -1:1, N, N);
    U = A\fx;

    % add on first and last nodes, U_0 and U_{N+1}
    U = [U(1)/(1 + h); U; U(end)/(1 + h)];

    % find exact solution
    ux = u(x)';

    % record error
    E = [E; h, norm(U - ux, inf)];
    
    plot(x,ux, x,U);
    pause;
end

% show ratio of decrease in error
E(1:end-1,2)./E(2:end,2)
