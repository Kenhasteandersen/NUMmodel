%
% Create an animation of a field from a global run.
%
% In:
%  x,y: ranges typically from sim.x and sim.y
%  field: with three dimensions: x, y, and time
%
% Example:
%  B = calcIntegrateGlobal(sim, sim.B);
%  animateGlobal(sim, log10(B), vContourLevels=[0 2],...
%      sProjection='ortho',bSpin=true, color=[0 0 0], bColorbar=false, time=20);
%
function F = animateGlobal(sim,field, options)

arguments
    sim;
    field (:,:,:) double = 0;
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
    
    options.bIncludeWatercolumn logical = false; % Use "plotGlobalWatercolum"
    options.lat double = 60; % The latitude to select
    options.lon double = -10; % The longitude to select
    
    options.width double = -1;
    options.height double = -1;
end

n = length(sim.t);
F(n) = struct('cdata',[],'colormap',[]);

figure(1)

if options.width ~= -1
    pos = get(gcf,'position')
    set(gcf,'Units','pixels')
    pos(3) = options.width;
    pos(4) = options.height;
end

for iTime = 1:n
    clf
    
    if options.bIncludeWatercolumn
        time = iTime/n*sim.t(end);
        plotGlobalWatercolumn(sim, time, options.lat, options.lon)
    else
        c = panelGlobal(sim.x,sim.y, field(:,:,iTime), options.vContourLevels, ...
            sTitle=options.sTitle, sProjection=options.sProjection);
        c.Label.String  = options.sUnits;
        
        %caxis([0 options.limit]);
        % Set bckground color
        set(gcf,'color',options.color);
        set(gca,'color',options.color);
        gridm('off'); % Remove grid lines
        framem('off'); % Remove lines around map
        
        % Possible delete the colorbar:
        if ~options.bColorbar
            colorbar off
        end
        
        % Do spin:
        if options.bSpin
            setm(gca,'Origin',[20 -iTime/n*360 0])
        end
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
v.Quality=95;

open(v);
for i = 1:n
    writeVideo(v,F(i));
end
close(v)

