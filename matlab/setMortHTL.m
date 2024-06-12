function setMortHTL(mortHTL)

setHTL(0,0,false,false);
calllib(loadNUMmodelLibrary(), 'f_setmorthtl', ...
    double(mortHTL));

% Set on parallel cluster as well:
if exist("gcp")
    if ~isempty(gcp('nocreate'))
        spmd
            setHTL(0,0,false,false);
            calllib(loadNUMmodelLibrary(), 'f_setmorthtl', ...
                double(mortHTL));
        end
    end
end
