function [U, Ux, Uy] = Poisson2D(f, g, L, N)
    p = inputParser;
    p.addRequired('f', @Utils.isFunctionHandle);
    p.addRequired('g', @Utils.isFunctionHandle);
    p.addRequired('L', @(x) isnumeric(x) && x > 0);
    p.addRequired('N', @Utils.isInteger);
    p.parse(f, g, L, N);

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
    F = h^2*f(Ux, Uy);
    %arrayfun(@(k) h^2*f(x(iFun(k) + 1), y(jFun(k) + 1)), k);
    % add on boundary conditions to F
    % add bottom boundary
    kBottom = 1:N;
    F(kBottom) = F(kBottom) + g(x(2:N+1),0);
    % add top boundary
    kTop = (N^2-N+1):N^2;
    F(kTop) = F(kTop) + g(x(2:N+1), L);
    % add left boundary
    kLeft = 1:N:(N^2-N)+1;
    F(kLeft) = F(kLeft) + g(0, y(2:N+1));
    % add right boundary
    kRight = N:N:N^2;
    F(kRight) = F(kRight) + g(L, y(2:N+1));

    % build sparse matrix A
    % build main diagonal
    iMain = k';
    jMain = k';
    sMain = 4*ones(N^2,1);
    % build upper main diagonal
    iUpper = k';
    iUpper(N:N:N^2) = [];
    jUpper = k';
    jUpper(1:N:N^2) = [];
    sUpper = -ones(size(iUpper)); 
    % build lower main diagonal
    iLower = jUpper;
    jLower = iUpper;
    sLower = sUpper;
    % build off upper diagonal
    iU = (1:(N^2-N))';
    jU = ((N+1):N^2)';
    sU = -ones(size(iU));
    % build off lower diagonal
    iL = jU;
    jL = iU;
    sL = sU;
    
    iA = [iMain; iUpper; iLower; iU; iL];
    jA = [jMain; jUpper; jLower; jU; jL];
    sA = [sMain; sUpper; sLower; sU; sL];
    A = sparse(iA, jA, sA);

    U = A\F';

    % change U from 1D vector to 2D matrix
    %U = vec2mat(U, N);
end
