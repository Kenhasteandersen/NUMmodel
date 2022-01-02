%
% Performs the testfunction in sEvaluation
%
function bSuccess = testFunction(sEvaluation, sDescription)
if nargin()==1
    sDescription = sEvaluation;
end


fprintf('Testing %s ... ', sDescription);

try
    bSuccess = eval(sEvaluation);
catch ME
    fprintf('Error: %s\n', ME.identifier);
    bSuccess = false;
end

if bSuccess
    fprintf('OK.\n');
else
    fprintf(2,'Failed!\n');
end
