f = @(x,y) -1.25*exp(x + .5*y);
u = @(x,y) exp(x + .5*y);
fLaplacian = @(x, y) -1.5625*exp(x + .5*y);
Error = [];
for N = [10*2.^(0:5)-1]
    [U, Ux, Uy] = Poisson2D_9PointLaplacian(f, u, 1, N, fLaplacian);
    uExact = u(Ux, Uy);
%     USquare = vec2mat(U,N);
%     uExactSquare = vec2mat(uExact,N);
%     figure;
%     surf(USquare);
%     hold on
%     surf(uExactSquare);
%     pause
    Error = [Error; Ux(2) - Ux(1), norm(U - uExact', inf)];
    
end
hRatios = Error(1:end-1,1)./Error(2:end,1);
errorRatios = Error(1:end-1,2)./Error(2:end,2);
order = log(errorRatios)./log(hRatios);
table(hRatios, errorRatios, order)
