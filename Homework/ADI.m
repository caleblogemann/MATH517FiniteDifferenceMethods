function [u, Ux, Uy, k] = ADI(N, T, f, g)
    h = 1/(N+1);
    x = @(i) i*h;
    y = @(j) j*h;

    % find k such that k = O(h)
    Nt = ceil(T/h);
    k = T/Nt;

    t = @(n) n*k;

    % create functions to swap between column-wise ordering and i,j ordering
    %kFun = @(i, j) i + (j-1)*N;
    iFun = @(k) floor((k-1)/N) + 1;
    jFun = @(k) mod(k-1,N)+1;

    % permutation matrix to change from natural row wise ordering to natural
    % column wise ordering or vice versa
    % uCol = P*uRow or uRow = P*uCol
    i = 1:N^2;
    j = repmat(0:N:N^2-N, 1, N) + kron(1:N,ones(1, N));
    P = sparse(i, j, ones(1, N^2));

    % matrix to represent diffusion in one direction
    % D*uRow = Dx2*U and D*uCol = Dy2*U
    e = ones(N, 1);
    tridiagonal = 1/h^2*spdiags([e, -2*e, e], [-1, 0, 1], N, N);
    D = kron(speye(N), tridiagonal);

    % N^2 by N^2 identity matrix
    I = speye(N^2);

    % initiate matrix to store u at times 1:Nt
    u = zeros(N^2, Nt+1);
    % create initial conditions
    % vector of x-values for column-wise ordering k = 1:N^2
    Ux = x(iFun(1:N^2));
    % vector of y-values for column-wise ordering k = 1:N^2
    Uy = y(jFun(1:N^2));
    % find initial values in column-wise ordering
    u(:, 1) = f(Ux, Uy);

    % create index arrays for boundary
    % k-indices for bottom boundary in column-wise ordering or
    % left boundary in row-wise ordering
    % ie when either x or y is zero
    zeroIndices = 1:N:N^2;
    % k-indices for top boundary in column-wise ordering or
    % right boundary in row-wise ordering
    % ie when either x or y is one
    oneIndices = N:N:N^2;
    for n = 1:Nt
        % u starts in column-wise ordering
        % First stage
        % (I + k/2 Dy2)*u
        rhs = (I + k/2*D)*u(:, n);
        % add boundary conditions for y-direction at time t = tn
        bottomBoundary = k/(2*h^2)*g(t(n-1), x(1:N), 0);
        topBoundary = k/(2*h^2)*g(t(n-1), x(1:N), 1);
        rhs(zeroIndices) = rhs(zeroIndices) + bottomBoundary';
        rhs(oneIndices) = rhs(oneIndices) + topBoundary';
        % change to row-wise ordering
        rhs = P*rhs;
        % add boundary conditions for x-direction at time t = tn + k/2
        leftBoundary = k/(2*h^2)*g(t(n-1)+k/2, 0, y(1:N));
        rightBoundary = k/(2*h^2)*g(t(n-1)+k/2, 1, y(1:N));
        rhs(zeroIndices) = rhs(zeroIndices) + leftBoundary';
        rhs(oneIndices) = rhs(oneIndices) + rightBoundary';
        % solve for uStar which approximates u at t = tn + k/2
        % uStar is in row-wise ordering
        uStar = (I - k/2*D)\rhs;

        % second stage
        % rhs row-wise ordering
        rhs = (I + k/2*D)*uStar;
        % add boundary conditions for x-direction at time t = tn + k/2
        % left and right boundaries same as before
        rhs(zeroIndices) = rhs(zeroIndices) + leftBoundary';
        rhs(oneIndices) = rhs(oneIndices) + rightBoundary';
        % change to column-wise ordering
        rhs = P*rhs;
        % add boundary conditions for y-direction at time t = tn + k
        bottomBoundary = k/(2*h^2)*g(t(n-1)+k, x(1:N), 0);
        topBoundary = k/(2*h^2)*g(t(n-1)+k, x(1:N), 1);
        rhs(zeroIndices) = rhs(zeroIndices) + bottomBoundary';
        rhs(oneIndices) = rhs(oneIndices) + topBoundary';
        % solve for u at time t = tn + k
        % u is in column wise ordering
        u(:,n+1) = (I - k/2*D)\rhs;
    end
end
