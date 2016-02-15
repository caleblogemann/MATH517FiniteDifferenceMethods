f = @(x,y) -1.25*exp(x + .5*y);
u = @(x,y) exp(x + .5*y);
for N = 10*2.^(1:5)
    [U, x, y] = Poisson2D(N, f, f);
    uExact = ;
    norm(U - uExact);
end
