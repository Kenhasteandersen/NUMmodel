% xlabel_names=[{'Phyto_{0.6}'},{'U'},{'U-Phyto_{0.6}'},{'M'},{'POM'},{'\SigmaB'}];
% legend_namesAll={'Pico','Nano','Micro'};
% mTitle='Biomass of plankton';


function stackedBarplot_percentage(BpnmVec_all,legend_namesAll,xlabel_names,mTitle)
newcolors = {'#b25781', '#66b9bb','#b5fded',' #b5eded','#b5dded',...
    '#b5cded','#b5bded', '#d7bde2'};



x2= 1:size(BpnmVec_all,1);
rng('default') % for reproducibility
y2 = BpnmVec_all';
h2 =  bar(x2, y2,'stacked'); 

% Compute percentage
yp2 = y2./sum(y2) * 100; 

% % Compute bar segment centers
xbarCnt2 = vertcat(h2.XEndPoints); 
ybarTop2 = vertcat(h2.YEndPoints);
ybarCnt2 = ybarTop2 - y2/2; 

% Create text strings
txt2 = compose('%.1f%%',yp2);

% Add text 
% figure()
% clf
th2 = text(xbarCnt2(:), ybarCnt2(:), txt2(:), ...
    'HorizontalAlignment', 'center', ....
    'VerticalAlignment', 'middle', ...
    'Color', 'w',....
    'FontSize', 8);
% xlabels
legend(legend_namesAll,'Location','bestoutside');
title(mTitle)
ylabel('Biomass (\mugC/L)')
xticklabels(xlabel_names)
% axis tight
colororder(newcolors);