%
% Calculate the radius of a particle from its mass. The calculation assumes
% a spherical particle with constant density.
%
% In:
%  m: mass in microgram carbon
%
% Out:
%  r: radius in micrometer
%
function r = calcRadius(m)

arguments
    m double;
end

rho = 0.4*1e6*1e-12; % mug/cm3 (Andersen et al 2016; rho = m/V = 0.3*d^3/(4/3*pi*(d/2)^3) )
r = (3/(4*pi)*m/rho).^(1/3);
