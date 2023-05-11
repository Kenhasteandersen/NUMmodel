%
% Search for a parameter value in a namelist file
%
function [parval] = search_namelist(filename,namelist,parameter)

fulltext = fileread(filename);

namelist = lower( strcat('input_',namelist) );
namelist( regexp(namelist,' ')) = '_';


TextAsCells = strip(regexp(fulltext, '\n', 'split'))';
nml_start_idx = find(strcmp(TextAsCells, {['&',namelist]}'));

if length(nml_start_idx)<1
    disp('No namelist by that name. Please check again')
    return
end

if length(nml_start_idx)>1
    disp('There are several namelists like that in the file')
    return
end

nml_end_idx = find(strcmp(TextAsCells', '/'));
i=nml_end_idx-nml_start_idx;
j=i==min(i(i>0));

the_area = TextAsCells(nml_start_idx:nml_end_idx(j));
mask = ~cellfun(@isempty, strfind(the_area, parameter));
the_one_line = the_area(mask);

newChr = split(the_one_line);
idx=find(strcmp(newChr, '='));

parval = cellfun(@str2num,newChr((idx+1)));
end