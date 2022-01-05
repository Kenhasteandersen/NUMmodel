testFunction('testLoadFortran');

testSetup('setupGeneralistsOnly', 2399);
testSetup('setupGeneralistsDiatoms_simple', 3204);
testSetup('setupGeneralistsDiatoms', 6034);
testSetup('setupDiatoms_simpleOnly', 781);
testSetup('setupGeneralists_cspOnly', 1366);
testSetup('setupGeneric', 2321);

testFunction('testChemostat');
testFunction('testChemostatEuler');
testFunction('testGlobal');