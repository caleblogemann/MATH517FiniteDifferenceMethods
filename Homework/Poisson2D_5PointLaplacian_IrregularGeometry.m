function [U, Ux, Uy, A] = Poisson2D_5PointLaplacian_IrregularGeometry(f, L, N)
    p = inputParser;
    p.addRequired('f', @Utils.isFunctionHandle);
    p.addRequired('L', @(x) isnumeric(x) && x > 0);
    p.addRequired('N', @Utils.isOddInteger);
    p.parse(f, L, N);

    h = L/(N+1);
    x = 0:h:L;
    y = 0:h:L;
    m = ceil(N/2);

    % indices not in geometry
    indices = N*(floor((0:m^2-1)/m) + N-m) + mod(0:m^2-1,m) + 1;

    % create functions to swap between row-wise ordering and i,j ordering
    %kFun = @(i, j) i + (j-1)*N;
    iFun = @(k) mod(k-1,N)+1;
    jFun = @(k) floor((k-1)/N) + 1;

    k = 1:N^2;
    % create vectors for x and y positions for each entry in U
    Ux = x(iFun(k) + 1);
    Uy = y(jFun(k) + 1);

    % vector of forcing function values at each index k
    F = h^2*f(Ux, Uy);
    F(indices) = 0;
    % boundary conditions are zero

    % build sparse matrix A
    % A is block tridiagonal with symmetric upper and lower diagonals
    e = ones(N, 1);
    % block on main diagonal
    T = spdiags([-e, 4*e, -e], [-1, 0, 1], N, N);
    % shape of main diagonal and off diagonals
    I = eye(N);
    O = spdiags([e, e], [-1, 1], N, N);
    A = kron(I, T) + kron(O, -I);

    % because of irregular geometry need to zero out indices in top left corner
    A(:, indices) = 0;
    A(indices, :) = 0;
    for i=indices
        A(i,i) = 1;
    end

    U = A\F';

    % change U from 1D vector to 2D matrix
    %U = vec2mat(U, N);
end
