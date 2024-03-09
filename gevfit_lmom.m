function phat = gevfit_lmom(data) 
    
lmomn=lmom(data,3);
tau3=lmomn(3)/lmomn(2);

% GEV
C = 2/(tau3 + 3)-log(2)/log(3);
kappa=7.8590*C+2.944*C^2;
beta=lmomn(2)*kappa/((1-2^(-kappa))*gamma(1+kappa));
alpha=lmomn(1)-beta/kappa*(1-gamma(1+kappa));
kappa = -kappa;
phat=[kappa beta alpha];

end

function [L] = lmom(X,nL)

% lmom by Kobus N. Bekker, 14-09-2004
% Based on calculation of probability weighted moments and the coefficient
% of the shifted Legendre polynomial.

% Given nonnegative integer nL, compute the 
% nL l-moments for given data vector X. 
% Return the l-moments as vector L.

[rows cols] = size(X);
if cols == 1 X = X'; end
n = length(X);
X = sort(X);
b = zeros(1,nL-1);
l = zeros(1,nL-1);
b0 = mean(X);

for r = 1:nL-1
    Num = prod(repmat(r+1:n,r,1)-repmat([1:r]',1,n-r),1);
    Den = prod(repmat(n,1,r) - [1:r]);
    b(r) = 1/n * sum( Num/Den .* X(r+1:n) );
end

tB = [b0 b]';
B = tB(length(tB):-1:1);

for i = 1:nL-1
    Spc = zeros(length(B)-(i+1),1);
    Coeff = [Spc ; LegendreShiftPoly(i)];
    l(i) = sum((Coeff.*B),1);
end

L = [b0 l];

end

function pk = LegendreShiftPoly(n)
if n==0 
    pk = 1;
elseif n==1
    pk = [2 -1]';
else
    pkm2 = zeros(n+1,1);
    pkm2(n+1) = 1;
    pkm1 = zeros(n+1,1);
    pkm1(n+1) = -1;
    pkm1(n) = 2;
    for k=2:n
        pk = zeros(n+1,1);
        for e=n-k+1:n
            pk(e) = (4*k-2)*pkm1(e+1) + (1-2*k)*pkm1(e) + (1-k)*pkm2(e);
        end
        pk(n+1) = (1-2*k)*pkm1(n+1) + (1-k)*pkm2(n+1);
        pk = pk/k;
        if k<n
            pkm2 = pkm1;
            pkm1 = pk;
        end 
    end
end

end