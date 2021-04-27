function p = setupDiatomsOnly(n)

arguments
    n {mustBeInteger, mustBePositive} = 10;
end

loadNUMmodelLibrary();
calllib(loadNUMmodelLibrary(), 'f_setupdiatomsonly', int32(n) );

p.n = n+3;
p.nGroups = 1;

p.idxN = 1;
p.idxDOC = 2;
p.idxSi = 3;
p.idxB = 4; % We have three nutrient groups so biomass groups starts at index 4.


p.n = n+p.idxB-1;
p.nGroups = 1;
p.typeGroups = 1;
p.ixStart = p.idxB;
p.ixEnd = p.n;

p = getMass(p);

p.u0(1:3) = [150, 0, 150]; % Initial conditions (and deep layer concentrations)
p.u0(p.idxB:p.n) = 1;