testFunction('testLoadFortran');

testSetup('setupGeneralistsOnly', 2397);
testSetup('setupGeneralistsDiatoms_simple', 1251);
testSetup('setupGeneralistsDiatoms', 6033);
testSetup('setupDiatoms_simpleOnly', -1170);
testSetup('setupGeneralists_cspOnly', 1366);
testSetup('setupGeneric', 2321);

testFunction('testChemostat');
testFunction('testChemostatEuler');
testFunction('testGlobal');