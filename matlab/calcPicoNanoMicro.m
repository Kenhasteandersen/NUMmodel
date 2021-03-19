%
% Returns the biomass of pico, nano and micro plankton in mugC/liter.
% Only works with generalists.
% Example: call calcPicoNanoMicro(sim.B(end,:), sim.p.pGeneralists);
%
function Bpnm = calcPicoNanoMicro(B,p)

m0 = min(p.m);
m2 = 2.38761e-06;
m20 = 0.00238761;
m200 = 2.38761;


Bpico = calcBiomassRange(B, p.m(3:end), m0,m2);
Bnano = calcBiomassRange(B, p.m(3:end), m2,m20);
Bmicro = calcBiomassRange(B, p.m(3:end), m20,m200);

Bpnm = [Bpico, Bnano, Bmicro];
