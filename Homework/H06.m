H = [];
E = [];
for m = [11, 23, 47, 95, 191, 383, 767, 1535]
    [h, k, err] = heat_CN(m);
    H = [H; h];
    E = [E; err];
end
hRatios = H(1:end-1)./H(2:end);
errorRatios = E(1:end-1)./E(2:end);
order = log(errorRatios)./log(hRatios);
table(hRatios, errorRatios, order)
