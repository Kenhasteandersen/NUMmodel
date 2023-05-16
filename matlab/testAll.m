
loadNUMmodelLibrary

testSetup('setupGeneralistsOnly',  3814);
testSetup('setupGeneralistsDiatoms_simple', 5487);
testSetup('setupGeneralistsDiatoms', 8409);
testSetup('setupDiatoms_simpleOnly', 1284);
testSetup('setupDiatomsOnly', 4623);% ???
testSetup('setupGeneric', 4545);

testFunction('testChemostat',96966);
testFunction('testChemostatEuler',1);
testFunction('testChemostatSeasonal',NaN);
testFunction('testGlobal',847346);