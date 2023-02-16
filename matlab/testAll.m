
loadNUMmodelLibrary

testSetup('setupGeneralistsOnly', 3916);
testSetup('setupGeneralistsDiatoms_simple', 4351);
testSetup('setupGeneralistsDiatoms', 7596);
testSetup('setupDiatoms_simpleOnly', 411);
!testSetup('setupGeneralists_cspOnly', 1444);
testSetup('setupGeneric', 4281);

testFunction('testChemostat',4348916);
testFunction('testChemostatEuler',1);
testFunction('testChemostatSeasonal',NaN);
testFunction('testGlobal',891937);