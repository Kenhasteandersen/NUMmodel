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

jR = rates.jR(ix);
set(gca,'yscale','linear')
fillbetweenlines(m, 0*jR, jR, [0 0 1]);
hold on

ymax = max(rates.jR);

switch p.nameGroup{iGroup}
    case 'Generalists'

        % Find beta parameters from the input file:
        betaL = search_namelist('../input/input.h', p.nameGroup{iGroup}, 'bL');
        betaN = search_namelist('../input/input.h', p.nameGroup{iGroup}, 'bN');
        betaDOC = search_namelist('../input/input.h', p.nameGroup{iGroup}, 'bDOC');
        betaF = search_namelist('../input/input.h', p.nameGroup{iGroup}, 'bF');
        betaG = search_namelist('../input/input.h', p.nameGroup{iGroup}, 'bg');

        jR_L = betaL*rates.jLreal(ix);
        jR_N = betaN*rates.jN(ix);
        jR_DOC = betaDOC*rates.jDOC(ix);
        jR_F = betaF*rates.jFreal(ix);
        jR_g = betaG*rates.f(ix).*rates.jMax(ix); % because jR_g= betaG*jNet = betaG*f*jMax
%         jR_g_old = betaG*rates.jTot(ix); this is slightly smaller than
%         jR_g
        jTot = rates.jTot(ix);% rates.f(ix).*rates.jMax(ix);
        
        fillbetweenlines(m, jR, jR+jR_DOC, [0 1 0]);
        fillbetweenlines(m, jR+jR_DOC, jR+jR_DOC+jR_L, [0 0.8 0]);
        fillbetweenlines(m, jR+jR_DOC+jR_L, jR+jR_DOC+jR_L+jR_N, [0 0.6 0]);
        fillbetweenlines(m, jR+jR_DOC+jR_L+jR_N, jR+jR_DOC+jR_L+jR_N+jR_F, [0 0.4 0]);
        fillbetweenlines(m, jR+jR_DOC+jR_L+jR_N+jR_F, jR+jR_DOC+jR_L+jR_N+jR_F+jR_g, [0.6 0 0]);
        fillbetweenlines(m, jR+jR_DOC+jR_L+jR_N+jR_F+jR_g, jR+jR_DOC+jR_L+jR_N+jR_F+jR_g+jTot, [0.8 0 0] );

        ymax = max(jR+jR_DOC+jR_L+jR_N+jR_F+jR_g+jTot);
        legend({'Basal','DOC','Light','Nutrients','Feeding','Growth','Total growth'})
    case 'Diatoms'
        % Find beta parameters from the input file:
        betaL = search_namelist('../input/input.h', p.nameGroup{iGroup}, 'bL');
        betaN = search_namelist('../input/input.h', p.nameGroup{iGroup}, 'bN');
        betaDOC = search_namelist('../input/input.h', p.nameGroup{iGroup}, 'bDOC');
        betaSi = search_namelist('../input/input.h', p.nameGroup{iGroup}, 'bSi');
        betaG = search_namelist('../input/input.h', p.nameGroup{iGroup}, 'bg');

        jR_L = betaL*rates.jLreal(ix);
        jR_N = betaN*rates.jN(ix);
        jR_DOC = betaDOC*rates.jDOC(ix);
        jR_Si = betaSi*rates.jSi(ix);
%         jR_g = betaG*rates.jTot(ix);
        jR_g = betaG*rates.f(ix).*rates.jMax(ix);
        jTot = rates.jTot(ix);%rates.f(ix).*rates.jMax(ix);

        fillbetweenlines(m, jR, jR+jR_DOC, [0 1 0]);
        fillbetweenlines(m, jR+jR_DOC, jR+jR_DOC+jR_L, [0 0.8 0]);
        fillbetweenlines(m, jR+jR_DOC+jR_L, jR+jR_DOC+jR_L+jR_N, [0 0.6 0]);
        fillbetweenlines(m, jR+jR_DOC+jR_L+jR_N, jR+jR_DOC+jR_L+jR_N+jR_Si, [0 0.4 0]);
        fillbetweenlines(m, jR+jR_DOC+jR_L+jR_N+jR_Si, jR+jR_DOC+jR_L+jR_N+jR_Si+jR_g, [0.6 0 0]);
        fillbetweenlines(m, jR+jR_DOC+jR_L+jR_N+jR_Si+jR_g, jR+jR_DOC+jR_L+jR_N+jR_Si+jR_g+max(0,jTot), [0.8 0 0] );

        ymax = 1.05*max(jR+jR_DOC+jR_L+jR_N+jR_Si+jR_g+max(0,jTot));
        legend({'Basal','DOC','Light','Nutrients','Silicate','Growth','Total growth'})
end

title(append('Respiration of ',lower(p.nameGroup{iGroup})))
%semilogx(m, rates.jTot,'k-','linewidth',2);
set(gca,'xscale','log')

xlim([min(p.m(p.idxB:end)) max(m)])
ylim([0 ymax])

xlabel('Cell mass ({\mu}g_C)')
ylabel('Respiration (day^{-1})')

hold off
%------------------------------------------------------------------
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

    end
%     time = datestr(clock,'YYYY_mm_dd_HH_MM_SS');
% 
% saveas(gcf,['panelResp_',time,'.png'])
% p.rand
end
