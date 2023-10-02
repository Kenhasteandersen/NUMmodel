%
% Sensitivity analysis of a Water column simulation for specific
% coordinates.
% First extract NPP from satellite data by running
% EXTRACT_files_NPP_latlon.m
%


name_png = 'sensitivityNPP_';
directory = 'C:\Users\ampap\OneDrive - Danmarks Tekniske Universitet\Trophic efficiency\Compare Biomass\Comparisons with data\';
    
saved_png = append(directory,name_png,'.png');

mHTL = 1;
% mortHTL = .15;
MortHTLs = 0.1;
bHTLdecline = true;
bHTLquadratic = true;   
runningTime = 15*365;

% newSinkingPOM = logspace(log10(0.1), log10(100), 8);
newSinkingPOM = 0.78622;%logspace(log10(0.4), log10(0.8), 8);
newRemin2= logspace(log10(0.01), log10(0.8), 8);

% newSinkingPOM = linspace(0.5, 100, 10);
% MortHTLs = [0.05, 0.1 0.2 0.4];

Ntotal = zeros(length(newSinkingPOM), runningTime);
NPP = zeros(length(newRemin2), length(newSinkingPOM), runningTime);
NPP_annual = zeros(length(newSinkingPOM));
%%
%--------------------------------------------
%        Initialize coordinates
%--------......................--------
        lat_to_find = 0;
        lon_to_find = -172;
%
%--------------------------------------------
for iRemin2 = 1:length(newRemin2)
       % which parameters needs to be changed
        paramToReplace={'remin2';'remin2'};

        % which input list do they belong to?
        InputListName={'input_generalists';'input_diatoms'};

        % what are the new values?
        % remin2=0.5;
        remin2 = newRemin2(iRemin2);
        remin2d=remin2;

        % change to cell array
        combinedvalues=[remin2,remin2d];
        newvalue=cell(length(combinedvalues),1);
    for i=1:length(combinedvalues)
        newvalue{i}=[num2str(combinedvalues(i)),'d0'];
    end

    % run the script
    substituteInputParameters(paramToReplace,InputListName,newvalue) 
    
    for i = 1:length(newSinkingPOM)

        p = parametersWatercolumn( setupNUMmodel );
        % setHTL(mortHTL, mHTL, bHTLquadratic, bHTLdecline)
        setSinkingPOM(p, newSinkingPOM(i))
        p.tEnd = runningTime;
        sim = simulateWatercolumn(p,lat_to_find,lon_to_find);
        simWCSeason = calcFunctionsTh(sim);

        Ntotal(i,:) = simWCSeason.Ntot;
        NPP(iRemin2,i,:) = simWCSeason.ProdNet; % (mgC / m2 / day)
        % NPP_annual(i) = simWCSeason.ProdNetAnnual;

    end
end

%%

NPP_cell = {};
NPP_cell_month_mean = {};

for i = 1:length(newRemin2)

    for j = 1:length(newSinkingPOM)

        matr_NPP = reshape(NPP(i,j,end-359:end), 30, 12);% every column is a month (30, 146)
        NPP_cell{j} = matr_NPP;
        NPP_cell_month_mean{i,j} = mean(matr_NPP, 1);

    end
end



%%
% load('NPP_extracted.mat')
figure(1)
clf(1)
tiledlayout(1,2)

nexttile
for i = 1:length(newRemin2)

    for j = 1:length(newSinkingPOM)


        % plot(NPP_cell_month_mean{i, j}, 'o-r', 'LineWidth', 1)% mgC/m2/day
        plot(NPP_cell_month_mean{1, 1}(end-11:end), 'o-r', 'LineWidth', 1)
        hold on
        plot(NPP_extracted(1,:), 'ob--')
        plot(NPP_extracted(2,:), 'og--')
        plot(NPP_extracted(3,:), 'om--')

legend('NUM','Eppley model', 'Standard VGPM', 'CAFE','Location','best')
xlabel('Time (month)')
ylabel('NPP (mgC / m^2 /day)')
 mTitle = append('Lat: ', string(lat_to_find), ', Lon: ', string(lon_to_find));
 my_title = append('remin2 = ', string(newRemin2(i)), ' / sinking = ', string(newSinkingPOM(j)));
        title(mTitle)

    end
end

nexttile
geoscatter(lat_to_find, lon_to_find, 'filled')

save('NPP_cell_month_mean.mat');

% exportgraphics(gcf,['C:\Users\ampap\OneDrive - Danmarks Tekniske Universitet\...' ...
    % 'Trophic efficiency\Compare Biomass\Comparisons with data\NPPsensitivity1.png'])

exportgraphics(gcf,[saved_png])       


%%
   
saved_png2 = append(directory,name_png,'_details.png');
% uplim=max(max(NPP_cell_month_mean),max(NPP_extracted));
                
figure(2)
clf(2)
t=tiledlayout(2,4)
for i = 1:length(newRemin2)

    for j = 1:length(newSinkingPOM)

nexttile                              % without this
        % plot(NPP_cell_month_mean{i, j}(end-11:end), 'o-r', 'LineWidth', 1)% mgC/m2/day
        plot(NPP_cell_month_mean{i, j}, 'o-r', 'LineWidth', 1)% mgC/m2/day
        % plot(NPP_cell_month_mean{i, j}(end-11:end), 'o-r', 'LineWidth', 1)
        hold on
        plot(NPP_extracted(1,:), 'ob--')
        plot(NPP_extracted(2,:), 'og--')
        plot(NPP_extracted(3,:), 'om--')

% legend('NUM','Eppley model', 'Standard VGPM', 'CAFE','Location','best')
xlabel('Time (month)')
ylabel('NPP (mgC / m^2 /day)')
 mTitle = append('Lat: ', string(lat_to_find), ', Lon: ', string(lon_to_find));
 my_title = append('remin2 = ', string(newRemin2(i)), ' / sinking = ', string(newSinkingPOM(j)));
        title(my_title)

    end
end
title(t,mTitle, 'FontSize', 24)
lg  = legend('NUM','Eppley model', 'Standard VGPM', 'CAFE','NumColumns',2); 
lg.Layout.Tile = 'South'; % <-- Legend placement with tiled layout


exportgraphics(gcf,[saved_png2])       

%%  GROUP CONTRIBUTIONS
% depth_layer=1;
day = sim.p.tEnd - 200;
% Bnum=squeeze(mean(sim.B(:,depth_layer,:),1)); % average Biomass at specific depth layer

Bnum=squeeze(sim.B(day,depth_layer,:)); % average Biomass at specific depth layer
figure
PicoNanoMicroBarplots(sim,Bnum,depth_layer)
%%
exportgraphics(gcf,[append(directory,name_png,'_layer',num2str(depth_layer),'_barplots.png')])       






