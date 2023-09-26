%
% Returns a vector of strategies determined by the uptake rates
% Second return argument is colors corresponding to each strategy
%
function [strategy, col] = calcTrophicStrategy(p, rates, bDrawStrategies)

arguments
    p, rates struct
    bDrawStrategies = true; % Whether to draw strategies in the background
end

rhoCN = 5.68;

col = zeros(length(rates.jL),3);

for i = 1:length(rates.jL)
    strategy{i} = 'Unknown';

    if rates.jN(i)*rhoCN > rates.jL(i)
        strategy{i} = 'Light limited';
        col(i,:) = [0.749,1,0.749];
    else
        strategy{i} = 'Nutrient limited';
        col(i,:) = [0.749,0.749,1];
    end

    if rates.jDOC(i) > rates.jLreal(i)
        strategy{i} = 'Osmoheterotroph';
        col(i,:) = [0.9,0.71,0.71];
    end

    if (rates.jFreal(i)/rates.jL(i) > 0.25) || (rates.jF(i) > rates.jN(i)*rhoCN)
        strategy{i} = 'Mixotroph';
        col(i,:) = [1,0.749,0.749];

    end

    if (rates.jNloss(i) > 1e-5) && (rates.jN(i) < rates.jF(i)/rhoCN)
        strategy{i} = 'Heterotroph';
        col(i,:) = [1,0.6,0.6];
    end

end

if bDrawStrategies
    %
    % Background color depending on trophic strategies:
    %
    iGroup = 1; %select the group for which the background is drawn
    ix = (p.ixStart(iGroup):p.ixEnd(iGroup));
    m = p.m(ix);
    color = col(ix-(p.idxB-1),:);
    colori = color(1,:);
    Xlim = calcXlim(p);
    Xmin = Xlim(1);
    Xmax = Xlim(2);
    rectangle(Position=[Xmin,0.0001,(m(2)+m(1))/2-Xmin,500], FaceColor=colori, EdgeColor=colori);
    hold on
    stratn = strategy(ix-(p.idxB-1));
    captionedstrat = [stratn(1)];
    color2caption = [colori];
    for i = 2:length(ix)-1
        colori=color(i,:);
        step=(m(i)+m(i+1))/2-(m(i-1)+m(i))/2;
        rectangle(Position=[(m(i-1)+m(i))/2,0.0001,step,500], FaceColor=colori, EdgeColor=colori);

        %color for the legend, removing duplicates
        if ~ismember(stratn(i), captionedstrat)
            captionedstrat=[captionedstrat,stratn(i)];
            color2caption=cat(1,color2caption,colori);
        end
    end

    colori=color(length(ix),:);
    rectangle(Position=[(m(length(ix)-1)+m(length(ix)))/2,0.0001,Xmax-(m(length(ix)-1)+m(length(ix)))/2,500], FaceColor=colori, EdgeColor=colori);
    if ~ismember(stratn(length(ix)), captionedstrat)
        captionedstrat=[captionedstrat,stratn(length(ix))];
        color2caption=cat(1,color2caption,colori);
    end

    % Empty plots for setting legend
    dum=[];
    for i=1:length(captionedstrat)
        dummyplot=plot(NaN, NaN, 's', 'MarkerSize', 10, 'MarkerFaceColor', color2caption(i,:));
        dum=[dum,dummyplot];
    end

    set(gca,'XScale', 'log', "Layer", "top")
end



% strategy = rep('Unknown', p$n)
%   strategy[r$jN*p$rhoCN>r$jL] = "Light limited"
%   strategy[r$jL>=r$jN*p$rhoCN] = "Nutrient limited"
%   strategy[r$jDOC > r$jLreal] = "Osmoheterotroph"
%   strategy[(r$jNloss>1e-5) & (r$JN<r$JF/p$rhoCN)] = "Heterotroph"
%   strategy[((r$jFreal/r$jL > 0.25) | (r$jF>r$jN*p$rhoCN) )& !strategy=="Heterotroph"] = "Mixotroph"
%   return(strategy)
% }
