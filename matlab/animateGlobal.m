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
    %options.limit double = max(field(:)); % max limit for the colorscale
    options.sTitle char = "";
    options.sUnits char = ""; % e.g.: "Concentration (\mug_C/l)";
    options.sFilename char = "Global";
    options.sProjection char = "fast"; % or use e.g. "ortho"
    options.bSpin logical = false; % whether to spin the globe (works best
                                   % with sProjection="ortho".
    options.time double = 10; % The duration of the animations (in seconds)
    options.vContourLevels = double([min(field(:)), max(field(:))]); % Passed to panelGlobal
    options.color double = [1 1 1]; % background color
    options.bColorbar logical = true; % Whether to mke the colorbar
end

n = size(field,3);
F(n) = struct('cdata',[],'colormap',[]);

figure(1)
parfor iTime = 1:n
    clf
    c = panelGlobal(x,y,field(:,:,iTime), options.vContourLevels, ...
        sTitle=options.sTitle, sProjection=options.sProjection);
    c.Label.String  = options.sUnits;
    %caxis([0 options.limit]);
    % Set bckground color
    set(gcf,'color',options.color);
    set(gca,'color',options.color);
    
    % Possible delete the colorbar:
    if ~options.bColorbar
        colorbar off
    end
    
    % Do spin:
    if options.bSpin
        setm(gca,'Origin',[20 -iTime/n*360 0])
    end

    % Get the frame.
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
v.FrameRate = n / options.time;

open(v);
for i = 1:n
    writeVideo(v,F(i));
end
close(v)

