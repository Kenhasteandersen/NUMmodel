%% Example with several parameters that needs to be changed
% which parameters needs to be changed
 paramToReplace={'remin2';'remin2'};

% which input list do they belong to?
InputListName={'input_generalists';'input_diatoms'};

% what are the new values?
remin2=0.01;
remin2d=remin2;

% change to cell array
combinedvalues=[remin2,remin2d];
newvalue=cell(length(combinedvalues),1);
for i=1:length(combinedvalues)
    newvalue{i}=[num2str(combinedvalues(i)),'d0'];
end

% run the script
substituteInputParameters(paramToReplace,InputListName,newvalue)

%% Example with 1 parameter that needs to be changed
% which parameter needs to be changed
 paramToReplace='palatability';
 % which input list does it belong to?
InputListName='input_POM';
% what is the new values?
palatability=0.1;
% change to cell array
newvalue=[num2str(palatability),'d0'];
substituteInputParameters(paramToReplace,InputListName,newvalue)