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

ix = (p.ixStart(iGroup):p.ixEnd(iGroup))-p.idxB+1;
m = p.m(ix+p.idxB-1);

jR = rates.jR;
set(gca,'yscale','linear')
fillbetweenlines(m, 0*jR, jR, [0 0 1]);
hold on

ymax = max(rates.jR);

switch p.nameGroup{iGroup}
    case 'Generalists'

        % Find beta parameters from the input file:
        betaL = search_namelist('../input/input.nlm', p.nameGroup{iGroup}, 'bL');
        betaN = search_namelist('../input/input.nlm', p.nameGroup{iGroup}, 'bN');
        betaDOC = search_namelist('../input/input.nlm', p.nameGroup{iGroup}, 'bDOC');
        betaF = search_namelist('../input/input.nlm', p.nameGroup{iGroup}, 'bF');
        betaG = search_namelist('../input/input.nlm', p.nameGroup{iGroup}, 'bg');

        jR_L = betaL*rates.jLreal;
        jR_N = betaN*rates.jN;
        jR_DOC = betaDOC*rates.jDOC;
        jR_F = betaF*rates.jFreal;
        jR_g = betaG*rates.jTot;

        fillbetweenlines(m, jR, jR+jR_DOC, [0 1 0]);
        fillbetweenlines(m, jR+jR_DOC, jR+jR_DOC+jR_L, [0 0.8 0]);
        fillbetweenlines(m, jR+jR_DOC+jR_L, jR+jR_DOC+jR_L+jR_N, [0 0.6 0]);
        fillbetweenlines(m, jR+jR_DOC+jR_L+jR_N, jR+jR_DOC+jR_L+jR_N+jR_F, [0 0.4 0]);
        fillbetweenlines(m, jR+jR_DOC+jR_L+jR_N+jR_F, jR+jR_DOC+jR_L+jR_N+jR_F+jR_g, [0.6 0 0]);

        ymax = max(jR+jR_DOC+jR_L+jR_N+jR_F+jR_g);
        legend({'Basal','DOC','Light','Nutrients','Feeding','Growth'})
    case 'Diatoms'
        % Find beta parameters from the input file:
        betaL = search_namelist('../input/input.nlm', p.nameGroup{iGroup}, 'bL');
        betaN = search_namelist('../input/input.nlm', p.nameGroup{iGroup}, 'bN');
        betaDOC = search_namelist('../input/input.nlm', p.nameGroup{iGroup}, 'bDOC');
        betaSi = search_namelist('../input/input.nlm', p.nameGroup{iGroup}, 'bSi');
        betaG = search_namelist('../input/input.nlm', p.nameGroup{iGroup}, 'bg');

        jR_L = betaL*rates.jLreal;
        jR_N = betaN*rates.jN;
        jR_DOC = betaDOC*rates.jDOC;
        jR_Si = betaSi*rates.jSi;
        jR_g = betaG*rates.jTot;

        fillbetweenlines(m, 0*jR, jR, [0 0 1]);
        fillbetweenlines(m, jR, jR+jR_DOC, [0 1 0]);
        fillbetweenlines(m, jR+jR_DOC, jR+jR_DOC+jR_L, [0 0.8 0]);
        fillbetweenlines(m, jR+jR_DOC+jR_L, jR+jR_DOC+jR_L+jR_N, [0 0.6 0]);
        fillbetweenlines(m, jR+jR_DOC+jR_L+jR_N, jR+jR_DOC+jR_L+jR_N+jR_Si, [0 0.4 0]);
        fillbetweenlines(m, jR+jR_DOC+jR_L+jR_N+jR_Si, jR+jR_DOC+jR_L+jR_N+jR_Si+jR_g, [0.6 0 0]);

        ymax = max(jR+jR_DOC+jR_L+jR_N+jR_Si+jR_g);
        legend({'Basal','DOC','Light','Nutrients','Silicate','Growth'})
end

semilogx(m, rates.jTot,'k-','linewidth',2);
set(gca,'xscale','log')

xlim([min(m) max(m)])
ylim([0 ymax])

xlabel('Cell mass ({\mu}g_C)')
ylabel('Respiration (day^{-1})')

hold off
%
% h = fillbetweenlines(x,y1,y2,color)
%
% Fill between y1 (lower) and y2 (upper).
%
function h = fillbetweenlines(x,y1,y2,color)

% Grey is default color if no color argument given:
if nargin==3
    color = 0.5*[1 1 1];
end

% possibly flip vectors:
x = reshape(x,1,length(x));
y1 = reshape(y1,1,length(y1));
y2 = reshape(y2,1,length(y2));

% If color is a scalar, make it into a grey-scale:
if length(color)==1
    color = color*[1 1 1];
end

% Find lower limit:
ymin = min( [ylim() y1] );
if strcmp(get(gca, 'yscale'),'log') && (ymin<=0)
    ymin = min(y1(y1>0));
    y1(y1<=0) = ymin;
end

x = [x x(end:-1:1)];
y = [y1 y2(end:-1:1)];

h=fill(x,y,color);
set(h,'edgecolor',color,'edgealpha',0);






