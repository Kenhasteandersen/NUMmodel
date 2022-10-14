%
% Setup matrix for sinking. Used internally by simulateWatercolumn and
% simulateGlobal.
%
function [Asink,p] = calcSinkingMatrix(p, sim, nGrid)
%
% Get sinking velocities from libNUMmodel:
%
p.velocity = 0*p.m;
p.velocity = calllib(loadNUMmodelLibrary(), 'f_getsinking', p.velocity);
p.idxSinking = find(p.velocity ~= 0); % Find indices of groups with sinking
%
% Set up the matrix:
%
Asink = zeros(p.n,nGrid,nGrid);
for i = 1:nGrid
    k = p.velocity/(sim.dznom(i)*p.dtTransport);
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