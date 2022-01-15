testFunction('testLoadFortran');

testSetup('setupGeneralistsOnly', 2535);
testSetup('setupGeneralistsDiatoms_simple', 3350);
testSetup('setupGeneralistsDiatoms', 6180);
testSetup('setupDiatoms_simpleOnly', 782);
testSetup('setupGeneralists_cspOnly', 1366);
testSetup('setupGeneric', 2471);

testFunction('testChemostat');
testFunction('testChemostatEuler');
testFunction('testGlobal');