%
% Performs the testfunction in sEvaluation
%
function bSuccess = testFunction(sEvaluation, value)
%if nargin()==1
sDescription = sEvaluation;
%end

fprintf('-----------------------------------\n');
fprintf('Testing %s ... ', sDescription);

try
    bSuccess = eval(sprintf('%s(%f)',sEvaluation,value));
catch ME
    fprintf('Error: %s\n', ME.identifier);
    bSuccess = false;
end

if bSuccess
    fprintf('OK.\n');
else
    fprintf(2,'Failed!\n');
end
