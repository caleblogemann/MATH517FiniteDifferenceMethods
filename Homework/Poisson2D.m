function [U] = Poisson2D(N, f, g)
    p = inputParser;
    p.addRequired('N', @Utils.isInteger);
    p.addRequired('f', @Utils.isFunctionHandle);
    p.addRequired('g', @Utils.isFunctionHandle);
    p.parse(N, f, g);

    h = 1/(N+1);
    x = 0:h:1;
    y = 0:h:1;

    % create functions to swap between row-wise ordering and i,j ordering
    kFun = @(i, j) i + (j-1)*N;
    iFun = @(k) mod(k,N);
    jFun = @(k) floor(k/N) + 1;

    % vector of forcing function values at each index k
    F = arrayfun(@(k) h^2*f(x(iFun(k) + 1), y(jFun(k) + 1)), 1:N^2);
    % add on boundary conditions to F
    G = zeros(N^2, 1);
    % add bottom boundary
    kBottom = find(jFun(1:N^2) == 1);
    G(kBottom) = g(x(kBottom));


    F = F + G;

    % build sparse matrix A
    T = spdiags([-1, 4, 1], -1:1, N);
    I = speye(N);
    A = sparse(N^2, N^2);
    for s = 1:N
        % put T as diagonal block
        A((s-1)*N+1:s*N,(s-1)*N+1:s*N) = T;
        if(s > 1)
            A((s-2)*N+1:(s-1)*N,(s-1)*N+1:s*N) = -I;
        end
        if(s < N)
            A(s*N+1:(s+1)*N,(s-1)*N+1:s*N) = -I;
        end
    end

    U = A\F;

    % change U from 1D vector to 2D matrix
    U = vec2mat(U, N);
end
