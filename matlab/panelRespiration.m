%
% Makes a panel of the various contributions to respiration for a group
% with beta-factors (generalists or diatoms).
%
% In:
%  p - the parameter structure
%  rates - rates
%  iGroup - the group number to be plotted
%
function panelRespiration(p,rates, iGroup)

% Find beta parameters from the input file:
betaL = search_namelist('../input/input.nlm', p.nameGroup{iGroup}, 'bL');
betaN = search_namelist('../input/input.nlm', p.nameGroup{iGroup}, 'bN');
betaDOC = search_namelist('../input/input.nlm', p.nameGroup{iGroup}, 'bDOC');
betaG = search_namelist('../input/input.nlm', p.nameGroup{iGroup}, 'bg');

jR_L = betaL*rates.jLreal;
jR_N = betaN*rates.jN;
jR_DOC = betaDOC*rates.jDOC;
jR_g = betaG*rates.jTot;
jR = rates.jR;

ix = (p.ixStart(iGroup):p.ixEnd(iGroup))-p.idxB+1;
m = p.m(ix+p.idxB-1);

hold off
semilogx(m, jR);
hold on
semilogx(m, jR+jR_g);
semilogx(m, jR+jR_g+jR_DOC);
semilogx(m, jR+jR_g+jR_DOC+jR_N);
semilogx(m, jR+jR_g+jR_DOC+jR_N+jR_L);

xlabel('Cell mass ({\mu}g_C)')
ylabel('Respiration (day^{-1})')






