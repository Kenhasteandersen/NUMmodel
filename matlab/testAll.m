
loadNUMmodelLibrary

testSetup('setupGeneralistsOnly',  3825);
testSetup('setupGeneralistsDiatoms_simple', 5513);
testSetup('setupGeneralistsDiatoms', 8426);
testSetup('setupDiatoms_simpleOnly', 1284);
testSetup('setupDiatomsOnly', 4623);% ???
testSetup('setupGeneric', 4563);

testFunction('testChemostat',157117);
testFunction('testChemostatEuler',1);
testFunction('testChemostatSeasonal',NaN);
testFunction('testGlobal',876467);