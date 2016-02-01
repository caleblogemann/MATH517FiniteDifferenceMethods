%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Problem 9
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% solving BVP -u'' + u = f(x) for f(x) = sin(4pi x) on 0 \le x \le 1
% with periodic BCs u(0) = u(1) and u'(0) = u'(1)

L = 1;
f = @(x) sin(4*pi*x);

% exact solution
u = @(x) 1/(16*pi^2 + 1)*sin(4*pi*x);

% create vector for storing errors at respective h values
E = [];
% N is number of interior point
for N=[10, 20, 40, 80, 160, 320, 640, 1280]
    % find spacing
    h = L/(N+1);

    % create x grid and function values on grid
    % x_0 to x_{N+1}
    x = 0:h:L;
    % leave out x_{N+1}
    fx = f(x(1:end-1))';

    % create matrix A that solves the finite difference problem
    e = ones(N+1, 1);
    A = 1/(12*h^2) * spdiags([-16*e, e, e, -16*e, (30 + 12*h^2)*e, -16*e, e, e, -16*e], [-N, -N+1, -2:2, N-1, N], N+1, N+1);
    U = A\fx;

    % add last node at x_{N+1}
    % U_{N+1} = U_{0}
    U = [U; U(1)];

    % compute exact solution
    ux = u(x)';

    % record error
    E = [E; h, norm(U - ux, inf)];

    % plot approximate and exact solution
    plot(x,ux, x,U);
    %pause;
end

% show ratio of decrease in error
hRatios = E(1:end-1,1)./E(2:end,1);
errorRatios = E(1:end-1,2)./E(2:end,2);
order = log(errorRatios)./log(hRatios);
table(hRatios, errorRatios, order)
