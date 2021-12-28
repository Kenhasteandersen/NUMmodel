%
% Performs the testfunction of the setup in sSetup
%
function bSuccess = testSetup(sSetup, value, sDescription)
arguments
    sSetup string;
    value double = 0;
    sDescription string = sSetup;
end

fprintf('Testing %s ... ', sDescription);

try
    p = eval( sSetup );
    u = p.u0;
    u(2) = 0.1;
    dudt = 0*u';
    [~, dudt] = calllib(loadNUMmodelLibrary(), 'f_calcderivatives', ...
        u, 60, 15, 0.1, dudt);
        
    if round(sum(dudt)*1000) == value
        bSuccess = true;
    else
        bSuccess = false;
        fprintf('sum(dudt) = %4.0f\n.',round(sum(dudt)*1000));
    end
catch ME
    fprintf('Error: %s\n', ME.identifier);
    bSuccess = false;
end

if bSuccess
    fprintf('OK.\n');
else
    fprintf('Failed!\n');
end
