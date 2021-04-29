function p = parametersAddgroup(typeGroup, p, n)

if ~isfield(p,'nGroups')
    p.nGroups = 1;
else
    p.nGroups = p.nGroups+1;
end
p.typeGroups(p.nGroups) = typeGroup;
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
ix = p.ixStart(p.nGroups):p.ixEnd(p.nGroups);

p.n = p.n + n;
