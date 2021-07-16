function p = parametersAddgroup(typeGroup, p, n)

if ~isfield(p,'nGroups')
    p.nGroups = 1;
    p.idxB = p.n+1;
else
    p.nGroups = p.nGroups+1;
end
p.typeGroups(p.nGroups) = typeGroup;

switch p.typeGroups(p.nGroups)
    case 1
        p.nameGroup{p.nGroups} = 'Generalists';
    case 2
        p.nameGroup{p.nGroups} = 'Generalists CSP';
    case 3
        p.nameGroup{p.nGroups} = 'Diatoms';
    case 4
        p.nameGroup{p.nGroups} = 'Diatoms simple';
    case 10
        p.nameGroup{p.nGroups} = 'Copepod';
end

%
% Set up indexing into u:
%
if (p.nGroups==1)
    iStart = p.idxB;
else
    iStart = p.ixEnd(p.nGroups-1)+1;
end
p.ixStart(p.nGroups) = iStart;
p.ixEnd(p.nGroups) = p.ixStart(p.nGroups)+n-1;
%ix = p.ixStart(p.nGroups):p.ixEnd(p.nGroups);

p.n = p.n + n;
