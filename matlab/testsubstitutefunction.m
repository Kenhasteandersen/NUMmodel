%% Example: change several parameters and write back to input.yaml

S = inputRead();

S.diatoms_simple.remin2 = 0.9;
S.generalists_simple.reminF = 0.8;

inputWrite(S);

%% Example: change a single parameter

S = inputRead();
S.diatoms_simple.remin2 = 0.9;
inputWrite(S);
