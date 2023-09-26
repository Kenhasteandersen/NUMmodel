%
% Makes a panel of the various contributions to exudation for a group
%
% In:
%  p - the parameter structure
%  rates - rates
%  iGroup - the group number to be plotted
%
function panelExudation(p,rates, iGroup)

ix = (p.ixStart(iGroup):p.ixEnd(iGroup))-p.idxB+1;
m = p.m(ix+p.idxB-1);

jLossPassive =(1-rates.f(ix)).* rates.jLossPassive(ix);
set(gca,'yscale','linear')
fillbetweenlines(m, 0*jLossPassive, jLossPassive, [0 0 1]);
hold on

ymax = max(rates.jR);

switch p.nameGroup{iGroup}
    case 'Generalists'

        eL = search_namelist('../input/input.h', p.nameGroup{iGroup}, 'epsilonL');
        eF = search_namelist('../input/input.h', p.nameGroup{iGroup}, 'epsilonF');

        jCloss_L=(1-eL)/eL*rates.jLreal(ix);
        jCloss_F=(1-eF)/eF*rates.jFreal(ix);

        jTot = rates.jTot(ix);%rates.f(ix).*rates.jMax(ix);
        
        fillbetweenlines(m, jLossPassive, jLossPassive+jCloss_L, [0 1 0]);
        fillbetweenlines(m, jLossPassive+jCloss_L, jLossPassive+jCloss_L+jCloss_F, [0 0.8 0]);
        fillbetweenlines(m, jLossPassive+jCloss_L+jCloss_F, jLossPassive+jCloss_L+jCloss_F+jTot, [0.8 0 0] )

        ymax = max(jLossPassive+jCloss_L+jCloss_F+jTot);
        legend({'Passive losses','Photouptake losses','Feeding losses','Total growth'})
    case 'Diatoms'
        
        eL = search_namelist('../input/input.h', p.nameGroup{iGroup}, 'epsilonL');
        
        jCloss_L=(1-eL)/eL*rates.jLreal(ix);

        jTot = rates.jTot(ix);%rates.f(ix).*rates.jMax(ix);
        
        fillbetweenlines(m, jLossPassive, jLossPassive+jCloss_L, [0 1 0]);
        fillbetweenlines(m, jLossPassive+jCloss_L, jLossPassive+jCloss_L+jTot, [0.8 0 0] )

        ymax = max(jLossPassive+jCloss_L+jTot);
        legend({'Passive losses','Photouptake losses','Total growth'})
end

title(append('Exudation of ',lower(p.nameGroup{iGroup})))
%semilogx(m, rates.jTot,'k-','linewidth',2);
set(gca,'xscale','log')

xlim([min(m) max(m)])
ylim([0 ymax])

xlabel('Cell mass ({\mu}g_C)')
ylabel('Respiration (day^{-1})')

hold off
%------------------------------------------------------------------
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

    end
%     time = datestr(clock,'YYYY_mm_dd_HH_MM_SS');
% 
% saveas(gcf,['panelResp_',time,'.png'])
% p.rand
end