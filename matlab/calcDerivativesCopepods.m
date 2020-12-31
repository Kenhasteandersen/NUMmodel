function rates = calcDerivativesCopepods(ixStart, ixEnd, u, rates, p)
ix = ixStart:ixEnd;
m = p.m(ix);
%
% Growth and reproduction:
%
nu = p.pCopepods.epsilonF.*rates.JF(ix) - p.pCopepods.Jresp;
g = max(zeros(1, length(ix)), nu)./m;
rates.mortStarve(ix) = -min(zeros(1,length(ix)), nu)./m;
b = p.pCopepods.epsilonR * g(end); % Birth rate
%
% Mortality:
%
rates.mort(ix) = rates.mortpred(ix) + p.mortHTLm(ix) + rates.mortStarve(ix);
%
% Assemble derivative:
% 
gamma = (g-rates.mort(ix)) ./ (1 - p.pCopepods.mZ.^(1-rates.mort(ix)./g));

rates.dudt(ixStart) = b*u(ixEnd) + (g(1)-gamma(1)-rates.mort(ixStart))*u(ixStart);
rates.dudt(ix(2:end-1)) = gamma(1:end-2).*u(ixStart:(ixEnd-2)) + ...
    (g(2:end-1)-gamma(2:end-1)-rates.mort(ix(2:end-1))).*u(ix(2:end-1));
rates.dudt(ixEnd) = gamma(end-1)*u(ixEnd-1) - rates.mort(ixEnd)*u(ixEnd);



rates.Jtot(ix) = nu;