testFunction('testLoadFortran');

testSetup('setupGeneralistsOnly', 2405);
testSetup('setupGeneralistsDiatoms_simple', 3221);
testSetup('setupGeneralistsDiatoms', 6052);
testSetup('setupDiatoms_simpleOnly', 782);
testSetup('setupGeneralists_cspOnly', 1366);
testSetup('setupGeneric', 2269);

testFunction('testChemostat');
testFunction('testChemostatEuler');
testFunction('testGlobal');