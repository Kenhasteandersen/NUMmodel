function p = parametersAddgroup(typeGroup, p, n, mass)
p.nGroups = p.nGroups+1;
p.typeGroups(p.nGroups) = typeGroup;
%
% Set up indexing into u:
%
if (p.nGroups==1)
    iStart = 3;
else
    iStart = p.ixEnd(p.nGroups-1)+1;
end
p.ixStart(p.nGroups) = iStart;
p.ixEnd(p.nGroups) = p.ixStart(p.nGroups)+n-1;
ix = p.ixStart(p.nGroups):p.ixEnd(p.nGroups);
%
% Get the parameters for the specific group:
%
switch typeGroup
    case p.typeGeneralists % Generalists
        pGroup = parametersGeneralists(n);
        p.pGeneralists = pGroup;
    case p.typeCopepods % Copepods
        pGroup = parametersCopepods(mass, n);
        p.pCopepods = pGroup;
    otherwise
        disp('Group number ',typeGroup,' not implemented')
        keyboard
end
%
% Extract grid into main parameter structure:
%
p.m = [p.m, pGroup.m];
p.mLower = [p.mLower, pGroup.mLower];
p.mDelta = [p.mDelta, pGroup.mDelta];
%
% Extract the feeding parameters into the main parameter structure
%

p.beta(ix) = pGroup.beta;
p.sigma(ix) = pGroup.sigma;

p.AF(ix) = pGroup.AF;
p.JFmax(ix) = pGroup.JFmax;
p.epsilonF(ix) = pGroup.epsilonF;
p.Jresp(ix) = pGroup.Jresp;

