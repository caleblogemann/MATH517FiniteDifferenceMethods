function [U, Ux, Uy] = Poisson2D_9PointLaplacian(f, g, L, N, fLaplacian)
    p = inputParser;
    p.addRequired('f', @Utils.isFunctionHandle);
    p.addRequired('g', @Utils.isFunctionHandle);
    p.addRequired('L', @(x) isnumeric(x) && x > 0);
    p.addRequired('N', @Utils.isInteger);
    p.addRequired('fLaplacian', @Utils.isFunctionHandle);
    p.parse(f, g, L, N, fLaplacian);

    h = L/(N+1);
    x = 0:h:L;
    y = 0:h:L;

    % create functions to swap between row-wise ordering and i,j ordering
    %kFun = @(i, j) i + (j-1)*N;
    iFun = @(k) mod(k-1,N)+1;
    jFun = @(k) floor((k-1)/N) + 1;

    k = 1:N^2;
    % create vectors for x and y positions for each entry in U
    Ux = x(iFun(k) + 1);
    Uy = y(jFun(k) + 1);

    % vector of forcing function values at each index k
    F = 6*h^2*(f(Ux, Uy)+h^2/12*fLaplacian(Ux, Uy));
    % add on boundary conditions to F
    % add bottom boundary
    kBottom = 1:N;
    F(kBottom) = F(kBottom) + 4*g(x(2:N+1),0) + g(x(1:N),0) + g(x(3:N+2),0);
    % add top boundary
    kTop = (N^2-N+1):N^2;
    F(kTop) = F(kTop) + 4*g(x(2:N+1), L) + g(x(1:N), L) + g(x(3:N+2),L);
    % add left boundary
    kLeft = 1:N:(N^2-N)+1;
    F(kLeft) = F(kLeft) + 4*g(0, y(2:N+1)) + g(0, y(1:N)) + g(0, y(3:N+2));
    % add right boundary
    kRight = N:N:N^2;
    F(kRight) = F(kRight) + 4*g(L, y(2:N+1)) + g(L, y(1:N)) + g(L, y(3:N+2));
    % each corner got a double boundary condition
    F(1) = F(1) - g(0, 0);
    F(N) = F(N) - g(L, 0);
    F(N^2 - N + 1) = F(N^2 - N + 1) - g(0, L);
    F(N^2) = F(N^2) - g(L, L);

    % build sparse matrix A
    % A is block tridiagonal with symmetric upper and lower diagonals
    e = ones(N, 1);
    % block on main diagonal
    T = spdiags([-4*e, 20*e, -4*e], [-1, 0, 1], N, N);
    % block for upper and lower diagonals
    S = spdiags([-1*e, -4*e, -1*e], [-1, 0, 1], N, N);
    % shape of main diagonals and off diagonals
    I = eye(N);
    O = spdiags([e, e], [-1, 1], N, N);
    A = kron(I, T) + kron(O, S);

    U = A\F';

    % change U from 1D vector to 2D matrix
    %U = vec2mat(U, N);
end
