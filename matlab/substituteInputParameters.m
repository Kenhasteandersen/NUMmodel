function substituteInputParameters(paramToReplace,InputListName,newvalue)
%this function creates a new input file with changed parameters based on
%the original input file in the input folder. The original input file is
%renamed from input.h to input_orig.h while the new input file is named
%input.h
%The function takes three inputs:
%
%  paramToReplace is the name of the parameter that will be replaced. This
%       can a string or cell array ex:
%       paramToReplace={'remin2';'reminF'}; or paramToReplace='remin2';
%
%  InputListName is the name of the input list the parameter belongs to. 
%       This can a string or cell array ex:
%       InputListName={'input_diatoms_simple';'input_generalists_simple'};
%       or InputListName='input_diatoms_simple';
%
%  newvalue is the new value of the parameter given as a string. This
%       can a string or cell array ex:
%       newvalue={'0.2d0','0.09d0'}; newvalue='0.2d0';

paramToReplace=cellstr(paramToReplace);
InputListName=cellstr(InputListName);
newvalue=cellstr(newvalue);

%% define paths and load lists
inputFileName='../input/input.h';
inputFileTmpName='../input/input_orig.h';
inputHeaderName='../input/input_parameter_headers.mat';
load(inputHeaderName,'thelists')
load(inputHeaderName,'thestrings')
%% move the file to a secure copy
status=movefile(inputFileName,inputFileTmpName);
if ~status
    error(['the input file ',inputFileName,' could not be renamed as ',inputFileTmpName])
end
% read in the input file
S = readlines(inputFileTmpName);
% find the lines where the different input lists begin
headerLineNr=nan(1,length(thelists));
for i=1:length(thelists)
headerLineNr(i) = find(contains(S,thestrings{i}));
end


for i=1:length(paramToReplace)
    whichline = find(contains(S,paramToReplace{i}));
    %which input file does it belong to
    whichInputList = strcmp(thelists,InputListName{i});
    theline=whichline(find(whichline>headerLineNr(whichInputList),1));
    thisline=S{theline};
    idx = find(thisline == '=');
    idx_comment = find(thisline == '!');
    idx_comment = min([idx_comment-1,length(thisline)]);
    value=strtrim(thisline(idx+1:idx_comment));
    S{theline} = regexprep(S{theline}, value, newvalue{i});
end

[fid, msg] = fopen(inputFileName, 'w');
if fid < 1
error('could not write output file because "%s"', msg);
end
fwrite(fid, strjoin(S, '\n'));
fclose(fid);
end