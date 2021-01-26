function rates = calcDerivativesCopepods(ixStart, ixEnd, u, rates, p)
ix = ixStart:ixEnd;
%
% Growth and reproduction:
%
nu = p.epsilonF*rates.JF(ix) - p.Jresp(ix);
g = max(zeros(1, length(ix)), nu)./p.m;
rates.mortStarve(ix) = min(zeros(1,length(ix)), nu)./p.m;
b = p.epsilonR * g(end); % Birth rate
%
% Mortality:
% p.mortHTL(ix)=p.mortHTL(ix).*u(ix);
%
rates.mort(ix) = rates.mortpred(ix) + p.mortHTL(ix) + rates.mortStart(ix);
%
% Assemble derivative:
% 
gamma = (g-rates.mort(ix)) ./ (1 - p.z(ix)).^(1-rates.mort(ix)./g);

rates.dudt(ixStart) = b*u(ixEnd) + (g(1)-gamma(1)-rates.mort(ixStart))*u(ixStart);
rates.dudt(ix(2:end-1)) = gamma(1:end-2)*u(ixStart:(ixEnd-2)) + ...
    (g(ix(2:end-1))-gamma(2:end-1)-rates.mort(ix(2:end-1))).*u(ix(2:end-1));
rates.dudt(ixEnd) = gamma(end-1)*u(ixEnd-1) - rates.mort(ixEnd)*u(ixEnd);