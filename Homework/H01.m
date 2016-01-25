%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Problem 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create test function
u = @(x) exp(x);
% Create exact second derivative
u2d = @(x) exp(x);

% randomly generate h
sz = [500,1];
H = 100*rand(sz);
h1 = H.*rand(sz);
h2 = H.*rand(sz);
h3 = H.*rand(sz);

% create x values from h values
x1 = -1*h1;
x2 = 0*H;
x3 = h2;
x4 = h2 + h3;

% compute weights
w1 = (2*(2*h2 + h3))./(h1.*(h1 + h2).*(h1 + h2 + h3));
w2 = (2*h1 - 4*h2 - 2*h3)./(h1.*h2.^2 + h1.*h2.*h3);
w3 = (2*(-1*h1 + h2 + h3))./(h2.*(h1 + h2).*h3);
w4 = (2*(h1 - h2))./(h3.*(h2 + h3).*(h1 + h2 + h3));

% approximate second derivative
u2da = w1.*u(x1) + w2.*u(x2) + w3.*u(x3) + w4.*u(x4);

% compute error
E = abs(u2d(x2) - u2da);
loglog(H, E, 'ko');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Problem 3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% do a least squares fit on the error in terms of H
A = [ones(sz), log(H)];
b = log(E);

% least squares fit Ax = b
lsf = (A'*A)\(A'*b);
lsf2 = A\b;
K = lsf(1);
p = lsf(2);
hold on
x = min(H):.01:max(H);
y = exp(K)*x.^p;
loglog(x,y);

