function test = testLayout(parentAxes)
%TESTLAYOUT Summary of this function goes here
%   Detailed explanation goes here
    test = tiledlayout(parentAxes, 3,1); % 2 lignes, 2 colonnes
    title(test, 'Mon Tiled Layout Complet');

    % Ajouter des sous-graphiques dans chaque tile
    ax1 = nexttile(test);
    plot(ax1, rand(1, 10));
    title(ax1, 'Graphique 1');

    ax2 = nexttile(test);
    bar(ax2, randi([1, 10], 1, 5));
    title(ax2, 'Graphique 2');

    ax3 = nexttile(test);
    pie(ax3, [1 2 3 4]);
    title(ax3, 'Graphique 3');
end

