%
% Returns the biomass of pico, nano and micro plankton in mugC/liter.
% Only works with generalists.
% Example: call calcPicoNanoMicro(sim.B(end,:), sim.p.pGeneralists);
%
function Bpnm = calcPicoNanoMicro(p,B)
% Calc mass from radius
rho = 0.4*1e6*1e-12; % mug/cm3 (Andersen et al 2016
mass = @(r) 4*pi/3*rho*r^3;
% Masses for pico, nano, and micro plankton
m0 = 0;
m2 = mass(1);
m20 = mass(10);
m200 = mass(100);

if size(B,1) ~= 1
    B = reshape(B,[1,length(B)]);
end

Bpico = sum( B.*calcMassRangeFraction(p, m0,m2) );
Bnano = sum( B.*calcMassRangeFraction(p, m2,m20) );
Bmicro = sum( B.*calcMassRangeFraction(p, m20,m200) );

Bpnm = [Bpico, Bnano, Bmicro];
