function [alphaF,palatability]=getPalatability
% parameters that define the vulnerability
gamma=0.6505;
v0_coefficient=1.7198;

alphaF=logspace(-4,0,20);
palatability=10^(v0_coefficient).*alphaF.^(gamma);
palatability=palatability./mean(palatability);
end