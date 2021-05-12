%
% Calculate the average rate of "j" across pico-nano-micro size groups
%
function jPNM = calcPicoNanoMicroRate(m,j)

m0 = min(m);
% upper pico-nano-micro sizes:
m2 = 2.38761e-06;
m20 = 0.00238761;
m200 = 2.38761;

% Average over pico-sizes
[~,f] = calcBiomassRange(j,m,m0,m2);
jPNM(1) = sum(f.*j)/sum(f);

% Average over nano-sizes
[~,f] = calcBiomassRange(j,m,m2,m20);
jPNM(2) = sum(f.*j)/sum(f);

% Average over micro-sizes
[~,f] = calcBiomassRange(j,m,m20,m200);
jPNM(3) = sum(f.*j)/sum(f);


