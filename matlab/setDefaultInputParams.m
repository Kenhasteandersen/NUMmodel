% sets default parameters to input.h
% to ensure we don't use a different version, e.g. after a sensitivity test
% the parameters are: remin, remin2, fracHTL_to_N, fracHTL_to_POM

function setDefaultInputParams 

paramToReplace={'remin2';'remin2';'remin';'fracHTL_to_N';'fracHTL_to_POM'};

% which input list do they belong to?
InputListName={'input_generalists';'input_diatoms';'input_POM';'input_general';'input_general'};

% what are the new values?
remin2=0.01;
remin2d=remin2;
remin=0.05;
fracHTL_to_N=0.5;
fracHTL_to_POM=0.5;

% change to cell array
combinedvalues=[remin2,remin2d,remin,fracHTL_to_N,fracHTL_to_POM];
newvalue=cell(length(combinedvalues),1);
for i=1:length(combinedvalues)
    newvalue{i}=[num2str(combinedvalues(i)),'d0'];
end

% run the script
substituteInputParameters(paramToReplace,InputListName,newvalue)

