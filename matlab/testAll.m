
loadNUMmodelLibrary

testSetup('setupGeneralistsSimpleOnly',  4170);
testSetup('setupGeneralistsOnly',  3727);
%testSetup('setupGeneralistsDiatoms_simple', 5487);
testSetup('setupGeneralistsDiatoms', 11380);
%testSetup('setupDiatoms_simpleOnly', 1284);
testSetup('setupDiatomsOnly', 7349);
%testSetup('setupGeneric', 4547);
testSetup('setupNUMmodel', 1630);

testFunction('testChemostat',363300);
testFunction('testChemostatEuler',1031.7);
%testFunction('testChemostatSeasonal',NaN);
testFunction('testWatercolumn',140121);
testFunction('testGlobal',1671930);