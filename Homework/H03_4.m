% Problem 6
f = @(x,y) ones(size(x));
for N = [19, 39]
    [U, Ux, Uy, A] = Poisson2D_5PointLaplacian_IrregularGeometry(f, 1, N);
    USquare = vec2mat(U,N);
    figure;
    surf(vec2mat(Ux, N), vec2mat(Uy, N), USquare);
    figure
    pcolor(USquare);
end

% Problem 7
f = @(x,y) 2*exp(-(10*x - 5).^2 - (10*y - 5).^2);
for N = [19, 39]
    [U, Ux, Uy, A] = Poisson2D_5PointLaplacian_IrregularGeometry(f, 1, N);
    USquare = vec2mat(U,N);
    figure;
    surf(vec2mat(Ux, N), vec2mat(Uy, N), USquare);
    figure
    pcolor(USquare);
end

% Problem 8 and 9
countsA = [];
countsB = [];
N = 10*2.^(0:5) - 1;
for iN = N
    [U, Ux, Uy, A] = Poisson2D_5PointLaplacian_IrregularGeometry(f, 1, iN);
    R = chol(A);
    % spy plots
    if(iN == 19 || iN == 39)
        figure
        spy(R);
        title(strcat('N= ',num2str(iN)));
    end
    countsA = [countsA; nnz(R)];
    P = symrcm(A);
    B = A(P, P);
    R = chol(B);
    if(iN == 19 || iN == 39)
        figure
        spy(R);
        title(strcat('Reverse Cuthill-Mckee, N= ',num2str(iN)));
    end
    countsB = [countsB; nnz(R)];
end
table(N', countsA)
table(N', countsB)
