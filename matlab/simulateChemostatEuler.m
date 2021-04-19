%
% Simulate the chemostat with Euler integration in Fortran library
% In:
%  p - parameter object (including chemostat parameters from
%      parametersChemostat). Not used if running from the fortran library. 
%  L - Light
% Out:
%  sim - simulation object
%
function sim = simulateChemostatEuler(p, L)
%
% Simulate:
%
u = p.u0;
dudt = 0*u;

u = calllib(loadNUMmodelLibrary(), 'f_simulatechemostateuler', ...
    int32(length(u)), u, L, p.u0(1), p.d, p.tEnd, 0.01);
%
% Assemble result:
%
sim.t = p.tEnd;
sim.N = u(1);
sim.DOC = u(2);
sim.B = u(3:end);
sim.p = p;
sim.rates = getRates(u(end,:), L);
for iGroup = 1:p.nGroups
    sim.Bgroup(:,iGroup) = sum( u(p.ixStart(iGroup):p.ixEnd(iGroup)));
end

end