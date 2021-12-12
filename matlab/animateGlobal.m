%
% Create an animation of a field from a global run. 
% Currently plots the plankton biomass
% and the nutrients.
%
% In: 
%  x,y: ranges typically from sim.x and sim.y
%  field: with three dimensions: x, y, and time
%
%
function F = animateGlobal(x,y,field, options)

arguments
    x,y (:,:) double;
    field (:,:,:) double;
    options.limit double = max(field(:)); % max limit for the colorscale
    options.sTitle char = "";
    options.sUnits char = ""; % e.g.: "Concentration (\mug C l^{-1})";
    options.sFilename char = "Global";
    options.sProjection char = "fast"; % or use e.g. "eckert4"
    options.bSpin logical = false; % whether to spin the globe (works best
                                   % with sProjection="ortho".
end

n = size(field,3);
F(n) = struct('cdata',[],'colormap',[]);

figure(1)
parfor iTime = 1:n
    clf
    c = panelGlobal(x,y,field(:,:,iTime), ...
        sTitle=options.sTitle, sProjection=options.sProjection);
    c.Label.String  = options.sUnits;
    caxis([0 options.limit]);
    
    if options.bSpin
        setm(gca,'Origin',[20 iTime/n*360 0])
    end

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

