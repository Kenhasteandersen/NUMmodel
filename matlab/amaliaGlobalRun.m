% runs a global simulation with the full NUMmodel setup    
    setDefaultInputParams 
    mHTL = 1;
    mortHTL = .15;
    bHTLdecline = true;
    bHTLquadratic = true;   
    newSinkingPOM = 1.54;

    p=setupNUMmodel();
    p = parametersGlobal(p); % Use standard low-res model
    p.tEnd = 5*365;
    setHTL(mortHTL, mHTL, bHTLquadratic, bHTLdecline)
    setSinkingPOM(p, newSinkingPOM)
    sim = simulateGlobal(p);
    sim.B(sim.B<0)=0; % Get rid of negative biomasses
    % checkConservation(sim);
