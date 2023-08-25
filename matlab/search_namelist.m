
% Search for a parameter value in a namelist file
%
function [parval] = search_namelist(filename,namelist,parameter)
load(fullfile('..','input','input_parameter_headers.mat'));

headernumber=find(strcmp(lower(namelist),erase(thelists,"input_"))==1);
if isempty(headernumber)
    disp('No namelist by that name. Please check again')
    return
elseif length(headernumber)>1
        disp('There are several namelists like that in the file')
        return
end
thislist='no list defined yet';

fid = fopen(filename,'r');
parval=[];
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
    % Check the parameter header is correct, check if the parameter exist
    %----------------------------------------------------------------------
    if strcmp(thislist,thelists(headernumber))
        %------------------------------------------------------------------
        % Remove leading and trailing white spaces
        %------------------------------------------------------------------
        line=strtrim(line);
        %------------------------------------------------------------------
        % If the line is not blank or a comment: extract the keyword 
        %------------------------------------------------------------------
        if ~isempty(line) && ~strcmp(line(1),'!')
            %--------------------------------------------------------------
            % Find the "=" and check if it is the parameter we search for
            %--------------------------------------------------------------
            idx = find(line == '=');
            keyword=strtrim(line(1:idx-1));
            %--------------------------------------------------------------
            % Find it is the right parameter, then extract the value
            %--------------------------------------------------------------
            if strcmp(keyword,parameter)
                idx_comment = find(line == '!');
                %----------------------------------------------------------
                % If there is no comment, use the length of the line 
                %----------------------------------------------------------
                idx_comment = min([idx_comment-1,length(line)]);
                value=strtrim(line(idx+1:idx_comment));
                value = strrep(value,'d','e');
                parval=str2double(value);
            end
        end
    end
end
if isempty(parval)
    disp(['Did not locate parameter: ',parameter])
end