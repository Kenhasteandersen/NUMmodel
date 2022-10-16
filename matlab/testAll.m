
loadNUMmodelLibrary

testSetup('setupGeneralistsOnly', 1174);
testSetup('setupGeneralistsDiatoms_simple', 1576);
testSetup('setupGeneralistsDiatoms', 4817);
testSetup('setupDiatoms_simpleOnly', 369);
testSetup('setupDiatomsOnly', 3619);% ???
testSetup('setupGeneralists_cspOnly', 1366);
testSetup('setupGeneric', 958);

testFunction('testChemostat',3888716);
testFunction('testChemostatEuler',1);
testFunction('testChemostatSeasonal',NaN);
testFunction('testGlobal',699480);