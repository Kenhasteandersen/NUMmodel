testFunction('testLoadFortran');

testSetup('setupGeneralistsOnly', 3767);
testSetup('setupGeneralistsDiatoms_simple', 4578);
testSetup('setupGeneralistsDiatoms', 7409);
testSetup('setupDiatoms_simpleOnly', 781);
testSetup('setupGeneralists_cspOnly', 1366);
testSetup('setupGeneric', 3837);

testFunction('testChemostat');
testFunction('testChemostatEuler');
testFunction('testGlobal');