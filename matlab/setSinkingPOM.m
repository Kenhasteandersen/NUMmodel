%
% Set the sinking velocity of POM
%
% Input:
%   p - the parameter struct
%   velocity - A vector with the length of the number POM size classes
%
function setSinkingPOM(p, velocity)

arguments
    p struct;
    velocity double;
end

if ~isfield( p, 'ixPOM' )
    error('There is no POM in this setup\n');
end

if length(velocity) ~= p.ixEnd(p.ixPOM)-p.ixStart(p.ixPOM)+1
    error('The length of the velocity vector is %i but should be %i\n', ...
        length(velocity), p.ixEnd(p.ixPOM)-p.ixStart(p.ixPOM)+1 );
end

vel = zeros(1,p.n); % All groups except POM has zero sinking
vel( p.ixStart(p.ixPOM):p.ixEnd(p.ixPOM) ) = velocity;

calllib(loadNUMmodelLibrary(), 'f_setsinking', vel );

% Set on parallel cluster as well:
if exist("gcp")
    if ~isempty(gcp('nocreate'))
        spmd
            calllib(loadNUMmodelLibrary(), 'f_setsinking', vel );
        end
    end
end