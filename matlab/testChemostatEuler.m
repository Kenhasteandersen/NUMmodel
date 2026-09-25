function bSuccess = testChemostatEuler(value)

p = setupNUMmodel;
p = parametersChemostat(p);
sim = simulateChemostatEuler(p, 100); % Same light level as testChemostat

bSuccess = true;

% Note: sim.B from simulateChemostatEuler is the final state, whereas
% simulateChemostat returns the whole time series, so the reference value
% here is not comparable with the one in testChemostat.
sumB = sum(sim.B(:));
if ~( sumB > 0.99*value && sumB < 1.01*value )
    bSuccess = false;
    fprintf(2,"sum(B) = %f\n",sumB);
end

fprintf("N balance = %+.2e 1/day. ", sim.Nbalance);
if abs(sim.Nbalance) > 1e-10
    bSuccess = false;
    fprintf(2,"N not conserved!\n");
end
