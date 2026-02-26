function S = inputRead(filename)

arguments
    filename char = '../input/input.yaml';
end
%--------------------------------------------------------------------------
% Create a structure for output
%--------------------------------------------------------------------------
S = struct();
thislist = '';
%--------------------------------------------------------------------------
% Open and read the YAML file, line by line
%
% Format:
%   section_name:          <- no leading space, ends with ':'
%     key: value  # comment  <- indented key-value pairs
%   # comment lines
%--------------------------------------------------------------------------
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
    % Key-value pair in the current section
    %----------------------------------------------------------------------
    if ~isempty(thislist)
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
        value = strtrim(trimline(idx+1:end));
        if isempty(value)
            continue
        end
        S.(thislist).(keyword) = str2double(value);
    end
end
fclose(fid);
