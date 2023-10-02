x = 1:5; 
rng('default') % for reproducibility
y = rand(4,5) * 10;  
h =  bar(x, y,'stacked'); 

% Compute percentage
yp = y./sum(y) * 100; 

% Compute bar segment centers
xbarCnt = vertcat(h.XEndPoints); 
ybarTop = vertcat(h.YEndPoints);
ybarCnt = ybarTop - y/2; 

% Create text strings
txt = compose('%.1f%%',yp);

% Add text 
th = text(xbarCnt(:), ybarCnt(:), txt(:), ...
    'HorizontalAlignment', 'center', ....
    'VerticalAlignment', 'middle', ...
    'Color', 'w',....
    'FontSize', 8);
%%
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
colororder(newcolors);

