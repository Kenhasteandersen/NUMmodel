
loadNUMmodelLibrary

testSetup('setupGeneralistsSimpleOnly',  4168);
testSetup('setupGeneralistsOnly',  3814);
testSetup('setupGeneralistsDiatoms_simple', 5487);
testSetup('setupGeneralistsDiatoms', 8409);
testSetup('setupDiatoms_simpleOnly', 1284);
testSetup('setupDiatomsOnly', 4623);% ???
testSetup('setupGeneric', 4545);
testSetup('setupNUMmodel', 1811);

testFunction('testChemostat',96966);
testFunction('testChemostatEuler',1);
%testFunction('testChemostatSeasonal',NaN);
testFunction('testWatercolumn',458403);
testFunction('testGlobal',10103073);