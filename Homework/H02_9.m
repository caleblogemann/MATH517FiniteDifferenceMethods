f = @(x,y) -1.25*exp(x + .5*y);
u = @(x,y) exp(x + .5*y);
Error = [];
for N = 10*2.^(0:5)
    [U, Ux, Uy] = Poisson2D(f, f, 1, N);
    uExact = u(Ux, Uy);
    Error = [Error, norm(U - uExact', 'inf')];
end
Error
