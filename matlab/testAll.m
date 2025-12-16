
loadNUMmodelLibrary

testSetup('setupGeneralistsSimpleOnly',  4170);
testSetup('setupGeneralistsOnly',  3727);
%testSetup('setupGeneralistsDiatoms_simple', 5487);
testSetup('setupGeneralistsDiatoms', 11354);
%testSetup('setupDiatoms_simpleOnly', 1284);
testSetup('setupDiatomsOnly', 7349);
testSetup('setupGeneric', 4547);
testSetup('setupNUMmodel', 1624);
testSetup('setupNUMmodelGelatinous',1393)

testFunction('testChemostat',26754);
testFunction('testChemostatEuler',1);
%testFunction('testChemostatSeasonal',NaN);
testFunction('testWatercolumn',189315);
testFunction('testGlobal',3071461);