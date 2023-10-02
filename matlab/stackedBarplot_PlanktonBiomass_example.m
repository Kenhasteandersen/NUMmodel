newcolors = {'#b25781', '#66b9bb','#b5fded',' #b5eded','#b5dded',...
    '#b5cded','#b5bded', '#d7bde2'};

 legend_namesVec=[{'Phyto_{0.6}'},{'U-Phyto_{0.6}'},{'M'},{'POM'}];


x2= 1:4;
rng('default') % for reproducibility
y2 = BpnmVec';
h2 =  bar(x2, y2,'stacked'); 

% Compute percentage
yp2 = y2./sum(y2) * 100; 


% Compute bar segment centers
xbarCnt2 = vertcat(h2.XEndPoints); 
ybarTop2 = vertcat(h2.YEndPoints);
ybarCnt2 = ybarTop2 - y2/2; 

% Create text strings
txt2 = compose('%.1f%%',yp2);

% Add text 
th2 = text(xbarCnt2(:), ybarCnt2(:), txt2(:), ...
    'HorizontalAlignment', 'center', ....
    'VerticalAlignment', 'middle', ...
    'Color', 'w',....
    'FontSize', 8);
legend({'Pico','Nano','Micro'},'Location','bestoutside');
title('Biomass of plankton')
ylabel('Biomass (\mugC/L)')
xticklabels(legend_namesVec)
colororder(newcolors);