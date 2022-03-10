loadNUMmodelLibrary

testSetup('setupGeneralistsOnly', 1215);
testSetup('setupGeneralistsDiatoms_simple', 1659);
testSetup('setupGeneralistsDiatoms', 4900);
testSetup('setupDiatoms_simpleOnly', 411);
testSetup('setupGeneralists_cspOnly', 1372);
testSetup('setupGeneric', 1049);

testFunction('testChemostat');
testFunction('testChemostatEuler');
testFunction('testGlobal');