%
% Create an animation of a field from a global run. 
% Currently plots the plankton biomass
% and the nutrients.
%
% In: 
%  field with three dimensions: x, y, and time
%
%
function F = animateGlobal(x,y,field, options)

arguments
    x,y (:,:) double;
    field (:,:,:) double;
    options.limit double = max(field(:));
    options.sTitle char = "";
    options.sUnits char = "Concentration (\mug N l^{-1})";
    options.sFilename char = "Global";
    options.sProjection char = "fast";
end

n = size(field,3);
F(n) = struct('cdata',[],'colormap',[]);

figure(1)
parfor iTime = 1:n
    clf
    c = panelGlobal(x,y,field(:,:,iTime),options.sTitle,options.sProjection);
    c.Label.String  = options.sUnits;
    caxis([0 options.limit]);

    drawnow
    F(iTime) = getframe(figure(1));
    disp(iTime);
end
%%
% Write the animation:
%
if (ismac || ispc)
    v = VideoWriter(options.sFilename, 'MPEG-4');
else
    v = VideoWriter(options.sFilename, 'Motion JPEG AVI');
end
v.FrameRate = n/5;

open(v);
for i = 1:n
    writeVideo(v,F(i));
end
close(v)

