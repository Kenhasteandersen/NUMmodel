%
% Setup with generalists and a number of copepods
%
function p = setupGeneric(mAdult, bParallel)

arguments
    mAdult (1,:) = [];
    bParallel = false;
end

loadNUMmodelLibrary();
calllib(loadNUMmodelLibrary(), 'f_setupgeneric', int32(length(mAdult)), mAdult );
if bParallel
    h = gcp('nocreate');
    poolsize = h.NumWorkers;
    parfor i=1:poolsize
        calllib(loadNUMmodelLibrary(), 'f_setupgeneric', int32(length(mAdult)), mAdult );
    end
end

p.idxN = 1;
p.idxDOC = 2;
p.idxB = 3; % We have two nutrient groups so biomass groups starts at index 3.

p.n = 2;
% Generalists:
p = parametersAddgroup(1,p,10);

for i = 1:length(mAdult)
    p = parametersAddgroup(10,p,10, mAdult(i));
end


%p.n = p.n+p.idxB-1;
%p.nGroups = 1+length(mAdult);
%p.typeGroups = ones(1,1+p.nGroups);
%p.typeGroups(2:end) = 10;
%p.ixStart = p.idxB;
%p.ixEnd = p.n;

p = getMass(p);

p.u0(1:2) = [150, 0]; % Initial conditions (and deep layer concentrations)
p.u0(p.idxB:p.n) = 1;