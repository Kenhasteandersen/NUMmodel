%
%
%
function p = parametersAddNutrients(p, nameShort, nameLong, nameUnits)

if ~isfield(p,'nNutrients')
    p.nNutrients=1;
else
    p.nNutrients = p.nNutrients + 1;
end

p.nameNutrientsShort{p.nNutrients} = nameShort;
p.nameNutrientsLong{p.nNutrients} = nameLong;
p.nameNutrientsUnits{p.nNutrients} = nameUnits;
p.(sprintf('idx%s',nameShort)) = p.nNutrients;
p.n = p.nNutrients;

if strcmp(p.nameNutrientsShort, 'N') 
    p.colNutrients{p.nNutrients} = [0 0 1];
end

if strcmp(p.nameNutrientsShort, 'DOC') 
    p.colNutrients{p.nNutrients} = [1 0 1];
end

if strcmp(p.nameNutrientsShort, 'Si') 
    p.colNutrients{p.nNutrients} = [0 1 1];
end

