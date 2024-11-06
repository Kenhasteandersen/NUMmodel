%
% Sets the HTL mortality to a specific vector (see setHTL for the
% "automatic" version).
%
% In:
%  mortHTL - the factor to multiply on pHTL
%  pHTL - the HTL selectivity
%  bQuadratic - whether the mortality is "quadratic" or not
%
function setMortHTL(mortHTL, pHTL, bQuadratic)

calllib(loadNUMmodelLibrary(), 'f_setmorthtl', ...
    double(mortHTL), double(pHTL), logical(bQuadratic));

% Set on parallel cluster as well:
if exist("gcp")
    if ~isempty(gcp('nocreate'))
        spmd
            calllib(loadNUMmodelLibrary(), 'f_setmorthtl', ...
                double(mortHTL), double(pHTL), logical(bQuadratic));
        end
    end
end
