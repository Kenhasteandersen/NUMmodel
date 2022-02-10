testFunction('testLoadFortran');

testSetup('setupGeneralistsOnly', 1001);
testSetup('setupGeneralistsDiatoms_simple', 1815);
testSetup('setupGeneralistsDiatoms', 4643);
testSetup('setupDiatoms_simpleOnly', 782);
testSetup('setupGeneralists_cspOnly', 1366);
testSetup('setupGeneric', 742);

testFunction('testChemostat');
testFunction('testChemostatEuler');
testFunction('testGlobal');