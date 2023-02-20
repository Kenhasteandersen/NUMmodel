function p = setupNutrients_N_DOC()

p = [];
p = parametersAddNutrients(p, 'N','Nitrogen','{\mug}N/l');
p = parametersAddNutrients(p, 'DOC','Dissolved organic matter','{\mug}C/l');
p.n=2;