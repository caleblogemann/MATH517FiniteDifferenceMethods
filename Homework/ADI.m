function [] = ADI(N, T, f, g)
    h = 1/(N+1);
    x = 0:h:1;
    y = 0:h:1;

    % find k such that k = O(h)
    Nt = ceil(T/h);
    k = T/Nt;

    % create functions to swap between row-wise ordering and i,j ordering
    %kFun = @(i, j) i + (j-1)*N;
    iFun = @(k) mod(k-1,N)+1;
    jFun = @(k) floor((k-1)/N) + 1;

    % permutation matrix to change from natural row wise ordering to natural
    % column wise ordering or vice versa
    % uCol = P*uRow or uRow = P*uCol
    i = 1:N^2;
    j = repmat(0:N:N^2-N, 1, N) + kron(1:N,ones(1, N));
    P = sparse(i, j, ones(1, N^2));

    % matrix to represent diffusion in one direction
    % D*uRow = Dx2*U and D*uCol = Dy2*U
    e = ones(N, 1);
    tridiagonal = spdiags([e, -2*e, e], [-1, 0, 1], N);
    D = kron(speye(N), tridiagonal);

    % N^2 by N^2 identity matrix
    I = speye(N^2);

    % create initial conditions

    for n = 1:Nt
        % u starts in column-wise ordering
        % First stage
        % (I + k/2 Dy2)*u
        rhs = (I + k/2*D)*u;
        % add boundary conditions for y-direction at time t = tn

        % change to row-wise ordering
        rhs = P*rhs;
        % add boundary conditions for x-direction at time t = tn + k/2

        % solve for uStar which approximates u at t = tn + k/2
        % uStar is in row-wise ordering
        uStar = (I - k/2*D)\rhs;

        % second stage
        % rhs row-wise ordering
        rhs = (I + k/2*D)*uStar
        % add boundary conditions for x-direction at time t = tn + k/2

        % change to column-wise ordering
        rhs = P*rhs;
        % add boundary conditions for y-direction at time t = tn + k

        % solve for u at time t = tn + k
        % u is in column wise ordering
        u = (I - k/2*D)\rhs;
    end
end
