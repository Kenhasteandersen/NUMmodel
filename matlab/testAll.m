loadNUMmodelLibrary

testSetup('setupGeneralistsOnly', 1677);
testSetup('setupGeneralistsDiatoms_simple', 2120);
testSetup('setupGeneralistsDiatoms', 5361);
testSetup('setupDiatoms_simpleOnly', 411);
testSetup('setupGeneralists_cspOnly', 1444);
testSetup('setupGeneric', 2042);

testFunction('testChemostat',3888716);
testFunction('testChemostatEuler',1);
testFunction('testChemostatSeasonal',NaN);
testFunction('testGlobal',699480);