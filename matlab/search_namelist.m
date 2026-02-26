
% Search for a parameter value in a YAML input file
%
% Inputs:
%   filename  - path to the YAML input file
%   namelist  - section name (e.g. 'general', 'generalists_simple')
%   parameter - parameter name to look up
%
function [parval] = search_namelist(filename, namelist, parameter)

thislist = '';
parval = [];

fid = fopen(filename, 'r');
while ~feof(fid)
    line = fgetl(fid);
    if ~ischar(line)
        continue
    end
    trimline = strtrim(line);
    %----------------------------------------------------------------------
    % Skip empty lines and comment lines
    %----------------------------------------------------------------------
    if isempty(trimline) || trimline(1) == '#'
        continue
    end
    %----------------------------------------------------------------------
    % Detect section headers: no leading whitespace, ends with ':'
    %----------------------------------------------------------------------
    if line(1) ~= ' '
        if trimline(end) == ':'
            thislist = trimline(1:end-1);
        end
        continue
    end
    %----------------------------------------------------------------------
    % Key-value pair: only process in the matching section
    %----------------------------------------------------------------------
    if strcmp(thislist, namelist)
        % Remove inline comments
        idx_comment = find(trimline == '#', 1);
        if ~isempty(idx_comment)
            trimline = strtrim(trimline(1:idx_comment-1));
        end
        % Find the first colon (key:value separator)
        idx = find(trimline == ':', 1);
        if isempty(idx)
            continue
        end
        keyword = strtrim(trimline(1:idx-1));
        if strcmp(keyword, parameter)
            value = strtrim(trimline(idx+1:end));
            parval = str2double(value);
        end
    end
end
fclose(fid);

if isempty(parval)
    disp(['Did not locate parameter: ', parameter])
end
