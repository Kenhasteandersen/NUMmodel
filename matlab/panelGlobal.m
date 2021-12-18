%
% Make a panel of a global field given in gridded coordinates.
%
% Input:
%  x,y: x and y coordinates (from sim.x and sim.y)
%  z: field to show (typically from sim)
%  vContourLevels: levels of contours, or just min, max values
% Optional:
%  sTitle: title on panel
%  sProjection: projection to use. Defaults to 'fast'. Other projections 
%               requires that the mapping toolbox is installed. 
%               Good projection is 'eckert4'.
%
% Out:
%  cbar: the color bar object
%
function cbar = panelGlobal(x,y,z, vContourLevels, options)

arguments
    x,y (:,1);
    z (:,:);
    vContourLevels = double([min(z(:)), max(z(:))]);
    options.sTitle string = '';
    options.sProjection string = 'fast';
end

% Adjust to global plot (close gap at lat 0)
z = [z;z(1,:)];
x = [x-x(1);360];

% Determine contour level if only min and max are given

if length(vContourLevels)==2
    vContourLevels = linspace(vContourLevels(1), vContourLevels(2),10);
end

z = double(squeeze( min(max(z,vContourLevels(1)),vContourLevels(end)))');

if (strcmp(options.sProjection,'fast'))
    contourf(x,y, z, vContourLevels, 'LineStyle','none');
    %contourf(x,y,squeeze(z)', vContourLevels, 'LineStyle','none');
    shading flat
    axis tight
else
    if exist('axesm') ~= 2
        error('Using projections requires that the mapping toolbox is installed.')
    end
    ax = axesm ( 'Origin',  [0 -90 0], 'MapProjection',options.sProjection, ...
        'Grid','on', 'Frame', 'on','ScaleFactor', 1, 'labelrotation',...
        'off', 'FLineWidth', 2);
    ax.XColor = 'white';
    ax.YColor = 'white';
    axis tight manual
    %plabel('PlabelLocation',20, 'PLabelMeridian', 91)
    contourfm(y,x,z, vContourLevels,'linestyle','none');
    %shading interp
    load coastlines
    h=patchm(coastlat,coastlon,0.4*[1 1 1]);
    set(h,'linestyle','none')
    %geoshow('landareas.shp', 'FaceColor', [0.8 0.8 0.8], 'EdgeColor',
    %'black'); % Slow
end

cbar = colorbar('eastoutside', 'FontSize',14);
cbar.Label.String  = '\mug C l^{-1}';
box off
title(options.sTitle)
caxis(vContourLevels([1,end]))


