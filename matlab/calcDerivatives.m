function rates = calcDerivatives(p,u, L)
%
% Split u
%
ix = 3:length(u);
B = max(zeros(1,length(u)-2), u(ix));
%
% Calc food available for all groups:
%
rates.F(ix) = (p.theta*B')';
rates.f(ix) = p.AF(ix).*rates.F(ix) ./ (p.AF(ix).*rates.F(ix) + p.JFmax(ix));
rates.JF(ix) = rates.f(ix) .* p.JFmax(ix);
%rates.JFreal = rates.JF;
%
% Calc resource uptake of unicellular groups:
%
for iGroup = 1:p.nGroups
    %switch p.typeGroups(iGroup)
    %    case 1 % Generalists
            rates = calcRatesGeneralists(p.ixStart(1),p.ixEnd(1), u, rates, p.pGeneralists, L);
    %end
end 
%
% Calc predation mortality for all groups:
%
rates.mortpred(3:length(u)) = ((p.theta') * (rates.JF(ix)./p.epsilonF(ix).*B./p.m(ix)./(rates.F(ix)+1e-100))')';
%
% Assemble derivatives
%
for iGroup = 1:p.nGroups
    switch(p.typeGroups(iGroup))
        case(1) % Generalists
            rates = calcDerivativesGeneralists(p.ixStart(iGroup), p.ixEnd(iGroup), u, rates, p);
        case(2) % Copepods
            rates = calcDerivativesCopepods(p.ixStart(iGroup), p.ixEnd(iGroup), u, rates, p);
    end
end

end
