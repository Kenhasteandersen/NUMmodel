%
% Make a panel of a global field given in gridded coordinates:
% If sForm = "fast" then it does a rough - but fast - map. Otherwise
% a nice map is drawn
function cbar = panelGlobal(x,y,z, sTitle, sProjection)
if (nargin==4)
    sProjection = 'fast';
end
% Adjust to global plot (close gap at lat 0)
z = [z;z(1,:)];
x = [x-x(1);360];

if (strcmp(sProjection,'fast'))
    surface(x,y,squeeze(z)');
    shading flat
    axis tight
else
    ax = axesm ( 'Origin',  [0 -90 0], 'MapProjection',sProjection, ...
        'Grid','on', 'Frame', 'on','ScaleFactor', 1, 'labelrotation',...
        'off', 'FLineWidth', 2);
    ax.XColor = 'white';
    ax.YColor = 'white';
    axis tight manual
    %plabel('PlabelLocation',20, 'PLabelMeridian', 91)
    surfacem(y,x ,squeeze(z)');
    %shading interp
    geoshow('landareas.shp', 'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'black');
end

cbar = colorbar('eastoutside', 'FontSize',14);
cbar.Label.String  = '\mug C l^{-1}';
box off
title(sTitle)


