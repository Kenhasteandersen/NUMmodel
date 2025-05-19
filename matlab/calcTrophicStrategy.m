function [strategy, col] = calcTrophicStrategy(p, rates, ax, bDrawStrategies)
arguments
    p, rates struct
    ax = []; % Optional target axis
    bDrawStrategies = true; % Whether to draw background rectangles
end



rhoCN = 5.68;
col = zeros(length(rates.jL), 3);
strategy = cell(length(rates.jL), 1);

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
    % --- Draw background color according to trophic strategy ---
    % No 'axes(ax)' here to avoid creating new figure windows
    hold(ax, 'on');

    iGroup = 1;
    ix = (p.ixStart(iGroup):p.ixEnd(iGroup));
    m = p.m(ix);
    color = col(ix - (p.idxB - 1), :);
    stratn = strategy(ix - (p.idxB - 1));

    Xlim = calcXlim(p);
    Xmin = Xlim(1);
    Xmax = Xlim(2);
    Ymin = 0.0001;
    Ymax = 500;

    colori = color(1,:);
    rectangle(ax, 'Position', [Xmin, Ymin, (m(2)+m(1))/2 - Xmin, Ymax - Ymin], ...
        'FaceColor', colori, 'EdgeColor', colori);

    captionedstrat = {stratn{1}};
    color2caption = colori;

    for i = 2:length(ix)-1
        colori = color(i,:);
        step = (m(i+1)+m(i))/2 - (m(i-1)+m(i))/2;
        rectangle(ax, 'Position', [(m(i-1)+m(i))/2, Ymin, step, Ymax - Ymin], ...
            'FaceColor', colori, 'EdgeColor', colori);

        % Add to legend if new strategy not already included
        if ~ismember(stratn{i}, captionedstrat)
            captionedstrat{end+1} = stratn{i};
            color2caption = [color2caption; colori];
        end
    end

    % Last rectangle
    colori = color(end,:);
    xStart = (m(end-1)+m(end))/2;
    rectangle(ax, 'Position', [xStart, Ymin, Xmax - xStart, Ymax - Ymin], ...
        'FaceColor', colori, 'EdgeColor', colori);

    if ~ismember(stratn{end}, captionedstrat)
        captionedstrat{end+1} = stratn{end};
        color2caption = [color2caption; colori];
    end

    % Dummy plots for legend
    dum = gobjects(1, length(captionedstrat));
    for i = 1:length(captionedstrat)
        dum(i) = plot(ax, NaN, NaN, 's', 'MarkerSize', 10, ...
            'MarkerFaceColor', color2caption(i,:), 'DisplayName', captionedstrat{i});
    end

    set(ax, 'XScale', 'log', "Layer", "top");
end
end
