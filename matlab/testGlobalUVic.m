%
% Test a global run on the UVic transport matrices (nTMmodel = 3).
%
% Checks the total biomass against a reference value and that nitrogen is
% conserved. UVic runs with a transport step of a third of a day and euler steps
% of a twelfth, so this also exercises simulateGlobal with a dtTransport and a dt
% that are not the 0.5 and 0.1 of the MITgcm matrices.
%
% In:
%  value - reference value of sum(sim.B) (default the value measured when this
%          test was written)
%
% Out:
%  bSuccess - whether the run reproduced the reference value and conserved N
%
function bSuccess = testGlobalUVic(value)

arguments
    value double = 3727901;
end

p = setupGeneralistsPOM(5,1); % Fast setup with POM
p = parametersGlobal(p, 3); % UVic transport matrices
p.tEnd = 30;
p.tSave = 10;

sim = simulateGlobal(p);

bSuccess = true;

sumB = sum(sim.B(~isnan(sim.B)));
if ~( sumB > 0.99*value && sumB < 1.01*value )
    bSuccess = false;
    fprintf(2,"sum(B) = %f\n",sumB);
end

[~, dNdt_per_N] = checkConservation(sim, false);
fprintf("relative N balance = %+.2e 1/yr. ", dNdt_per_N);
if abs(dNdt_per_N) > 1e-4
    bSuccess = false;
    fprintf(2,"N not conserved!\n");
end
