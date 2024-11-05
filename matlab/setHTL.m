
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
%  mHTL - the size where the HTL starts (usually set to the maximum size
%         divided by 500^1.5)
%  bQuadraticHTL - boolean which determines whether to use a fixed or
%                  "quadratic" mortality
%  bDecliningHTL - boolean which determines whether the HTL mortality
%                  declines with size as mass^-0.25.
%
% Out:
%  Nothing; the function only affects the fortran library.
%
function setHTL(mortHTL, mHTL, bQuadraticHTL, bDecliningHTL)

arguments
    mortHTL double {mustBeNonnegative} = 0.2;
    mHTL double = 1/500^1.5; % Suits simulations with only generalists
    bQuadraticHTL logical = false;
    bDecliningHTL logical = false;
end


calllib(loadNUMmodelLibrary(), 'f_sethtl', ...
    double(mHTL), double(mortHTL), logical(bQuadraticHTL), logical(bDecliningHTL) );

% Set on parallel cluster as well:
if exist("gcp")
    if ~isempty(gcp('nocreate'))
        spmd
            calllib(loadNUMmodelLibrary(), 'f_sethtl', ...
                double(mHTL), double(mortHTL), logical(bQuadraticHTL), logical(bDecliningHTL) );
        end
    end
end