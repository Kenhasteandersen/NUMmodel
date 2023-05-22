function S = inputRead(filename)

arguments
    filename char = '../input/input.h';
end
%--------------------------------------------------------------------------
% Creat a structure for output
%--------------------------------------------------------------------------
S = struct();
%--------------------------------------------------------------------------
% Define the parameter headers
%--------------------------------------------------------------------------
load(fullfile('..','input','input_parameter_headers.mat'));
% str_General='! GENERAL PARAMETERS';
% str_Generalists_simple='! GENERALISTS SIMPLE INPUT PARAMETERS';
% str_Generalists='! GENERALISTS INPUT PARAMETERS';
% str_Diatoms_simple='! DIATOMS SIMPLE INPUT PARAMETERS';
% str_Diatoms='! DIATOMS INPUT PARAMETERS';
% str_Copepods_active='! COPEPODS ACTIVE INPUT PARAMETERS';
% str_Copepods_passive='! COPEPODS PASSIVE INPUT PARAMETERS';
% str_POM='! PARTICULATE ORGANIC MATTER (POM) INPUT PARAMETERS';
% 
% 
% thelists={'input_general';'input_generalists_simple';'input_generalists';...
%     'input_diatoms_simple';'input_diatoms';'input_copepods_active';'input_copepods_passive';...
%     'input_POM'};
% 
% thestrings={str_General;str_Generalists_simple;str_Generalists;...
%     str_Diatoms_simple;str_Diatoms;str_Copepods_active;...
%     str_Copepods_passive;str_POM};

thislist='no list defined yet';
%--------------------------------------------------------------------------
% Open and read the text file, line by line
%--------------------------------------------------------------------------
fid = fopen(filename,'r');
while ~feof(fid)
    line = fgetl(fid);
    %----------------------------------------------------------------------
    % Check if this line is a parameter header. If, then update "thislist"
    %----------------------------------------------------------------------
    io=find(strcmp(line,thestrings)==1);
    if ~isempty(io)
        thislist=string(thelists(io,1));
    end
    %----------------------------------------------------------------------
    % Remove leading and trailing white spaces
    %----------------------------------------------------------------------
    line=strtrim(line);
    %----------------------------------------------------------------------
    % If the line is not blank or a comment: extract the keyword and value
    %----------------------------------------------------------------------
    if ~isempty(line) && ~strcmp(line(1),'!')
        %------------------------------------------------------------------
        % Find the "=" and "!" marking comments in the end of the line
        %------------------------------------------------------------------
        idx = find(line == '=');
        idx_comment = find(line == '!');
        %------------------------------------------------------------------
        % If there is no comment, use the length of the line instead
        %------------------------------------------------------------------
        idx_comment = min([idx_comment-1,length(line)]);
        %------------------------------------------------------------------
        % Seperate keyword from value (and convert to double)
        %------------------------------------------------------------------        
        keyword=strtrim(line(1:idx-1));
        value=strtrim(line(idx+1:idx_comment));
        value = strrep(value,'d','e');
        value=str2double(value);
        %------------------------------------------------------------------
        % Save to the current parameter list
        %------------------------------------------------------------------   
        S.(thislist).(keyword)=value;
    end
end
fclose(fid);

