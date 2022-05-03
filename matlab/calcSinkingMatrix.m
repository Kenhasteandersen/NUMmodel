%
% Setup matrix for sinking. Used internally by simulateWatercolumn and
% simulateGlobal.
%
function [Asink,p] = calcSinkingMatrix(p, sim, nGrid)

% Parameters for size dependent sinking speed
p.wRef=3.52;        % [m/d] Reference particle sinking velocity (Kriest & Oschlies (2008))
p.rRef=10;          % [µm]  Initial particle size
p.gamma=1.3;        % [-] Size dependent sinking velocity (Cram 2018)

% Parameter for depth dependent remineralisation
p.active=9*10^12;   % [molC/year] activation energy (John et al 2014) << find other [UNITS]

% Size dependent sinking velocity equation
p.velocity = 0*p.m;
p.velocity(13:22) = p.wRef*((calcRadius(p.m(13:22))/p.rRef).^p.gamma);
p.idxSinking = find(p.velocity ~= 0); % Find indices of groups with sinking

% Set up the matrix:
%
Asink = zeros(p.n,nGrid,nGrid);
for i = 1:nGrid
    k =p.velocity'/(sim.dznom(i)*p.dtTransport); % Forsøgte med en simpel ligning før vi kører på med den store.

    % ----- Den ligning vi tænker at prøve med
    %(p.velocity'*exp(-p.active*sim.z(i)*p.gamma))./(sim.dznom(i)*p.dtTransport); % Remineralisaton with depth affecting p.velocity
    % -----
    
    Asink(:,i,i) = 1+k;
    if (i ~= 1)
        Asink(:,i,i-1) = -k;
    end
end
%BottomBC: 
Asink(end) = 1;

%
% Invert matrix to make it ready for use:
%
for i = 1:p.n
    Asink(i,:,:) = inv(squeeze(Asink(i,:,:)));
end