function inputWrite(S, filename)

arguments
    S struct;
    filename char = '../input/input.yaml';
end
% inputWrite(S, FILENAME) writes a parameter structure S to a YAML file
% that can be read back by inputRead.m.
%
% S must follow the structure:
%
%     |-- section_a -- | -- var1
%     |                | -- var2
% S --|
%     |-- section_b -- | -- var1
%                      | -- var2
%
% Only scalar numeric values are supported (matching the NUM model input
% parameter format).

fid = fopen(filename, 'w');
sections = fieldnames(S);

for i = 1:length(sections)
    fprintf(fid, '%s:\n', sections{i});
    params = S.(sections{i});
    param_names = fieldnames(params);
    for j = 1:length(param_names)
        val = params.(param_names{j});
        fprintf(fid, '  %s: %g\n', param_names{j}, val);
    end
    fprintf(fid, '\n');
end

fclose(fid);
