function setMortHTL(mortHTL)

setHTL(0,0,false,false);
calllib(loadNUMmodelLibrary(), 'f_setmorthtl', ...
    double(mortHTL));
