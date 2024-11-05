%
% Calculate the average rate of "j" across pico-nano-micro size groups
%
function jPNM = calcPicoNanoMicroRate(p,j)

% Calc mass from radius
rho = 0.4*1e6*1e-12; % mug/cm3 (Andersen et al 2016
mass = @(r) 4*pi/3*rho*r^3;
% Masses for pico, nano, and micro plankton
m0 = 0;
m2 = mass(1);
m20 = mass(10);
m200 = mass(100);

% Average over pico-sizes
f = calcMassRangeFraction(p,m0,m2);
jPNM(1) = sum(f.*j)/sum(f);

% Average over nano-sizes
f = calcMassRangeFraction(p,m2,m20);
jPNM(2) = sum(f.*j)/sum(f);

% Average over micro-sizes
f = calcMassRangeFraction(p,m20,m200);
jPNM(3) = sum(f.*j)/sum(f);