function p = parametersAddgroup(typeGroup, p, n, mAdult)

arguments
    typeGroup ;
    p struct,
    n ;
    mAdult = 0;
end

if ~isfield(p,'nGroups')
    p.nGroups = 1;
    p.idxB = p.n+1;
else
    p.nGroups = p.nGroups+1;
end
p.typeGroups(p.nGroups) = typeGroup;

switch p.typeGroups(p.nGroups)
    case 1
        p.nameGroup{p.nGroups} = 'Generalists simple';
        p.colGroup{p.nGroups} = [0 0 0.75];
    %case 2
    %    p.nameGroup{p.nGroups} = 'Generalists CSP';
    %    p.colGroup{p.nGroups} = [0 0 0.75];
    case 3
        p.nameGroup{p.nGroups} = 'Diatoms';
        p.colGroup{p.nGroups} = [0 0.5 0];
    case 4
        p.nameGroup{p.nGroups} = 'Diatoms simple';
        p.colGroup{p.nGroups} = [0 0.5 0];
    case 5
        p.nameGroup{p.nGroups} = 'Generalists';
        p.colGroup{p.nGroups} = [0 0 0.75];
    case 10
        p.nameGroup{p.nGroups} = sprintf('Passive copepod %.1f {\\mu}g',mAdult);
        p.colGroup{p.nGroups} = [0.6 0.0 0];
    case 11
        p.nameGroup{p.nGroups} = sprintf('Active copepod %.1f {\\mu}g',mAdult);
        p.colGroup{p.nGroups} = [0.85 0.0 0];
    case 100
        p.nameGroup{p.nGroups} = 'POM';
        p.colGroup{p.nGroups} = [165 42 42]/256; % Brown
        p.ixPOM = p.nGroups;
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
