%
% Makes a panel of the various contributions to respiration for a group
% with beta-factors (generalists or diatoms).
%
% In:
%  p - the parameter structure
%  rates - rates
%  iGroup - the group number to be plotted
%  bScaled - flag to tell if the rates should be scaled by the total rates
%
function panelRespirationExudation(p,rates, iGroup, bScaled)

arguments
    p struct;
    rates;
    iGroup;
    bScaled logical = false;
end

ix = (p.ixStart(iGroup):p.ixEnd(iGroup))-p.idxB+1;
m = p.m(ix+p.idxB-1);

jR = rates.jR(ix);
jLossPassive =(1-rates.f(ix)).* rates.jLossPassive(ix);
set(gca,'yscale','linear')

ymax = max(rates.jR);

switch p.nameGroup{iGroup}
    case 'Generalists'

        % Find beta and epsilon parameters from the input file:
        betaL = search_namelist('../input/input.h', p.nameGroup{iGroup}, 'bL');
        betaN = search_namelist('../input/input.h', p.nameGroup{iGroup}, 'bN');
        betaDOC = search_namelist('../input/input.h', p.nameGroup{iGroup}, 'bDOC');
        betaF = search_namelist('../input/input.h', p.nameGroup{iGroup}, 'bF');
        betaG = search_namelist('../input/input.h', p.nameGroup{iGroup}, 'bg');
        eL = search_namelist('../input/input.h', p.nameGroup{iGroup}, 'epsilonL');
        eF = search_namelist('../input/input.h', p.nameGroup{iGroup}, 'epsilonF');

        jR_L = betaL*rates.jLreal(ix);
        jR_N = betaN*rates.jN(ix);
        jR_DOC = betaDOC*rates.jDOC(ix);
        jR_F = betaF*rates.jFreal(ix);
        jR_g = betaG*rates.f(ix).*rates.jMax(ix); % because jR_g= betaG*jNet = betaG*f*jMax
        jTot = rates.jTot(ix);
        jCloss_L=(1-eL)/eL*rates.jLreal(ix);
        jCloss_F=(1-eF)/eF*rates.jFreal(ix);

        if bScaled
            jTotal = jR+jR_DOC+jR_L+jR_N+jR_F+jR_g+jTot+jLossPassive+jCloss_L+jCloss_F;
            jR = jR ./ jTotal;
            jR_L = jR_L ./ jTotal;
            jR_N = jR_N ./ jTotal;
            jR_DOC = jR_DOC ./ jTotal;
            jR_F = jR_F ./ jTotal;
            jR_g = jR_g ./ jTotal;
            jTot = jTot ./ jTotal;
            jCloss_L = jCloss_L ./ jTotal;
            jCloss_F = jCloss_F ./ jTotal;
        end
        
        % Respirations:
        fillbetweenlines(m, 0*jR, jR, [0 0.2 0]);  % Basal
        hold on
        fillbetweenlines(m, jR, jR+jR_DOC, [0 1 0]); % DOC
        fillbetweenlines(m, jR+jR_DOC, jR+jR_DOC+jR_L, [0 0.8 0]); % Light
        fillbetweenlines(m, jR+jR_DOC+jR_L, jR+jR_DOC+jR_L+jR_N, [0 0.6 0]); % Nutrients
        fillbetweenlines(m, jR+jR_DOC+jR_L+jR_N, jR+jR_DOC+jR_L+jR_N+jR_F, [0 0.4 0]); % Feeding
        fillbetweenlines(m, jR+jR_DOC+jR_L+jR_N+jR_F, jR+jR_DOC+jR_L+jR_N+jR_F+jR_g, [0 0.3 0]); % Growth
        % Growth:
        fillbetweenlines(m, jR+jR_DOC+jR_L+jR_N+jR_F+jR_g, jR+jR_DOC+jR_L+jR_N+jR_F+jR_g+jTot, [0.8 0 0] ) % Growth
        % Exudations:
        fillbetweenlines(m, jR+jR_DOC+jR_L+jR_N+jR_F+jR_g+jTot, jR+jR_DOC+jR_L+jR_N+jR_F+jR_g+jTot+jLossPassive, [0 0 0.1]); % Passive
        fillbetweenlines(m, jR+jR_DOC+jR_L+jR_N+jR_F+jR_g+jTot+jLossPassive, jR+jR_DOC+jR_L+jR_N+jR_F+jR_g+jTot+jLossPassive+jCloss_L, [0 0 0.6]); % Light
        fillbetweenlines(m, jR+jR_DOC+jR_L+jR_N+jR_F+jR_g+jTot+jLossPassive+jCloss_L, jR+jR_DOC+jR_L+jR_N+jR_F+jR_g+jTot+jLossPassive+jCloss_L+jCloss_F, [0 0 0.4]); % Feeding

        ymax = max(jR+jR_DOC+jR_L+jR_N+jR_F+jR_g+jTot+jLossPassive+jCloss_L+jCloss_F);
        legend({'Basal','DOC','Light','Nutrients','Feeding','Growth','Total growth',...
            'Passive losses','Photouptake losses','Feedig losses'})
    case 'Diatoms'

        % Find beta and epsilon parameters from the input file:
        betaL = search_namelist('../input/input.h', p.nameGroup{iGroup}, 'bL');
        betaN = search_namelist('../input/input.h', p.nameGroup{iGroup}, 'bN');
        betaDOC = search_namelist('../input/input.h', p.nameGroup{iGroup}, 'bDOC');
        betaSi = search_namelist('../input/input.h', p.nameGroup{iGroup}, 'bSi');
        betaG = search_namelist('../input/input.h', p.nameGroup{iGroup}, 'bg');
        eL = search_namelist('../input/input.h', p.nameGroup{iGroup}, 'epsilonL');

        jR_L = betaL*rates.jLreal(ix);
        jR_N = betaN*rates.jN(ix);
        jR_DOC = betaDOC*rates.jDOC(ix);
        jR_Si = betaSi*rates.jSi(ix);
        %         jR_g = betaG*rates.jTot(ix);
        jR_g = betaG*rates.f(ix).*rates.jMax(ix);
        jTot = rates.jTot(ix);%rates.f(ix).*rates.jMax(ix);
        jCloss_L=(1-eL)/eL*rates.jLreal(ix);

        fillbetweenlines(m, jR, jR+jR_DOC, [0 1 0]);
        fillbetweenlines(m, jR+jR_DOC, jR+jR_DOC+jR_L, [0 0.8 0]); %
        fillbetweenlines(m, jR+jR_DOC+jR_L, jR+jR_DOC+jR_L+jR_N, [0 0.6 0]);
        fillbetweenlines(m, jR+jR_DOC+jR_L+jR_N, jR+jR_DOC+jR_L+jR_N+jR_Si, [0 0.4 0]); % Silcate uptake
        fillbetweenlines(m, jR+jR_DOC+jR_L+jR_N+jR_Si, jR+jR_DOC+jR_L+jR_N+jR_Si+jR_g, [0 0.3 0]);

        fillbetweenlines(m, jR+jR_DOC+jR_L+jR_N+jR_Si+jR_g, jR+jR_DOC+jR_L+jR_N+jR_Si+jR_g+jTot, [0.8 0 0] )

        fillbetweenlines(m, jR+jR_DOC+jR_L+jR_N+jR_Si+jR_g+jTot, jR+jR_DOC+jR_L+jR_N+jR_Si+jR_g+jTot+jLossPassive, [0 0 0.1]);
        fillbetweenlines(m, jR+jR_DOC+jR_L+jR_N+jR_Si+jR_g+jTot+jLossPassive, jR+jR_DOC+jR_L+jR_N+jR_Si+jR_g+jTot+jLossPassive+jCloss_L, [0 0 0.6]);

        ymax = max(jR+jR_DOC+jR_L+jR_N+jR_Si+jR_g+jTot+jLossPassive+jCloss_L);
        legend({'Basal','DOC','Light','Nutrients','Silicate','Growth','Total growth','Passive losses','Photouptake losses'})
end

title(append('Respiration and exudation of ',lower(p.nameGroup{iGroup})))
%semilogx(m, rates.jTot,'k-','linewidth',2);
set(gca,'xscale','log')

xlim([min(m) max(m)])
ylim([0 ymax])

xlabel('Cell mass ({\mu}g_C)')

if bScaled
    ylabel('Relative respiration or exudation')
else
    ylabel('Respiration or exudation (day^{-1})')
end

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