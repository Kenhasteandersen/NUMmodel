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