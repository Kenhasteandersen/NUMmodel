%
% Sets the higher trophic level mortality. Already done by default in the 
% setup function, so this function is only needed to override the defaults.
%
% HTL mortality is either a standard mortality that is fixed by mortHTL, or
% it is a "quadradratic" mortality that is proportional to the biomass. In
% the latter case morthHTL is the constant in front of the quadratic
% mortality, and not the mortality itself.
%
% In:
%  mortHTL - the HTL mortality (or the constant if bQuadratic=true)
%  bQuadraticHTL - boolean which determines whether to use a fixed or
%                  "quadratic" mortality
%
% Out:
%  Nothing; the function only affects the fortran library.
%
function parametersHTL(mortHTL, bQuadraticHTL)

arguments
    mortHTL double {mustBeNonnegative} = 0.2;
    bQuadraticHTL logical = false;
end

calllib(loadNUMmodelLibrary(), 'f_parametershtl', double(mortHTL), logical(bQuadraticHTL) );