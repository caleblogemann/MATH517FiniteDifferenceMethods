function [f] = uxFun(u, a, h)
    % f = 
    f = -a/(6*h)*([u(end-1); u(end); u(1:end-2)] - 6*[u(end); u(1:end-1)] + 3*u + 2*[u(2:end); u(1)]);
end
