%
% Returns the biomass of pico, nano and micro plankton in mugC/liter.
% Only works with generalists.
% Example: call calcPicoNanoMicro(sim.B(end,:), sim.p.pGeneralists);
%
function Bpnm = calcPicoNanoMicro(B,m)

m0 = min(m);
m2 = 2.38761e-06;
m20 = 0.00238761;
m200 = 2.38761;


Bpico = calcBiomassRange(B, m, m0,m2);
Bnano = calcBiomassRange(B, m, m2,m20);
Bmicro = calcBiomassRange(B, m, m20,m200);

Bpnm = [Bpico, Bnano, Bmicro];
