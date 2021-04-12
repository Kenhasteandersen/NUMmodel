%
% Create an animation of a global run. Currently plots the plankton biomass
% and the nutrients.
%
function F = animateGlobal(sim, sFilename, sProjection)

arguments
    sim struct;
    sFilename char = "Global";
    sProjection char = "fast";
end


n = length(sim.t);
F(n) = struct('cdata',[],'colormap',[]);

figure(1)
for iTime = 1:n
    clf
    tiledlayout(2,1)
    nexttile
    % Nutrients:
    c = panelGlobal(sim.x,sim.y,sim.N(:,:,1,iTime),'N',sProjection);
    caxis([0 350])
    c.Label.String  = 'Concentration [\mug N l^{-1}]';

    % Unicellular plankton
    nexttile
    panelGlobal(sim.x,sim.y,log10(sum(sim.B(:,:,1,findIxUnicellular(sim.p),iTime),4)),'Unicellular plankton (log10)',sProjection);
    caxis([1 3])
    ax2 = gca;
    ax2.NextPlot = 'replace';

    drawnow
    F(iTime) = getframe(figure(1));
    disp(iTime);
end
%%
% Write the animation:
%
if (ismac || ispc)
    v = VideoWriter(sFilename, 'MPEG-4');
else
    v = VideoWriter(sFilename, 'Motion JPEG AVI');
end
v.FrameRate = n/5;

open(v);
for i = 1:n
    writeVideo(v,F(i));
end
close(v)

