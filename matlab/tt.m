mAdult=[]
%p = setupGeneralistsOnly();

p = setupGeneralistsDiatoms_simple();
%   
% Set parameters:setupGeneralistsDiatoms_simple
%
    % n = 10;
    % nCopepods = 10;
    % nPOM = 10;
    % p = setupNUMmodel(mAdult, n,nCopepods,nPOM);
%


% p = setupNUMmodel()


p = parametersChemostat(p);
p.tEnd = 200;
p.d = 0.1;
%
% Set to "normal" HTL mortality if there are no copepods:
%
if isempty(mAdult)
    setHTL(0.1, 1/500^1.5, false, false);
else 
    setHTL(0.1, 1, true, true);
end
%
% Simulate
%
tic
sim = simulateChemostat(p, 100);
toc
%
% Plot
%
plotSimulation(sim);
checkConservation(sim);
