%
% Performs the testfunction of the setup in sSetup
%
function bSuccess = testSetup(sSetup, value, sDescription)
arguments
    sSetup string;
    value double = 0;
    sDescription string = sSetup;
end

fprintf('Testing %s ... \n', sDescription);

try  
    %
    % First check that the derivative function respects conservation
    % by making a short simulation
    %
    p = eval( sSetup ); % Initialize the setup
    p = parametersChemostat(p);
    p.tEnd = 1;
    sim = simulateChemostat(p);
    %plotSimulation(sim);
    fprintf('N balance: %e\n', checkBalanceDerivative(sim));
    %
    % Then check the value of the derivative:
    %
    p = eval( sSetup ); % Initialize the setup
    u = p.u0;
    u(2) = 0.1;
    dudt = 0*u';
    [~, dudt] = calllib(loadNUMmodelLibrary(), 'f_calcderivatives', ...
        u, 60, 15, 0.1, dudt);
        
    if round(sum(dudt)*1000) == value
        bSuccess = true;
    else
        bSuccess = false;
        fprintf(2,'sum(dudt) = %4.0f\n.',round(sum(dudt)*1000));
    end
catch ME
    fprintf(2,'Error: %s\n', ME.identifier);
    bSuccess = false;
end

if bSuccess
    fprintf('OK.\n');
else
    fprintf(2,'Failed!\n');
end
