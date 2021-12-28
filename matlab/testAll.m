testFunction('testLoadFortran');

testSetup('setupGeneralistsOnly', 1915);
testSetup('setupGeneralistsDiatoms_simple', 288);
testSetup('setupGeneralistsDiatoms', 5070);
testSetup('setupDiatoms_simpleOnly', -1652);
testSetup('setupGeneralists_cspOnly', 1366);
testSetup('setupGeneric', 2321);

testFunction('testChemostat');
testFunction('testChemostatEuler');
testFunction('testGlobal');