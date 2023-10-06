%
% Sensitivity analysis of a Water column simulation for specific
% coordinates.
% First extract NPP from satellite data by running
% EXTRACT_files_NPP_latlon.m
%
% load('NPP_extracted_seasonal.mat')
% lat_to_find=55; lon_to_find=-40;name_png='sensitivityNPP';noYears=20;

function Sensitivity_paramRemin2(lat_to_find,lon_to_find,name_png,NPP_extracted,noYears)
   
saved_png = append(name_png,'.png');

mHTL = 1;
mortHTL = .15;
bHTLdecline = true;
bHTLquadratic = true;   
runningTime = noYears*365;
newSinkingPOM = 0.78622;
newRemin2= logspace(log10(0.01), log10(0.5), 10);

Ntotal = zeros(length(newRemin2), runningTime);
NPP = zeros(length(newRemin2), runningTime);
NPP_annual = zeros(length(newRemin2));

for  iparam = 1:length(newRemin2)
       % which parameters needs to be changed
        paramToReplace={'remin2';'remin2'};

        % which input list do they belong to?
        InputListName={'input_generalists';'input_diatoms'};

        % what are the new values?
        remin2 = newRemin2(iparam);
        remin2d=remin2;

        % change to cell array
        combinedvalues=[remin2,remin2d];
        newvalue=cell(length(combinedvalues),1);
    for i=1:length(combinedvalues)
        newvalue{i}=[num2str(combinedvalues(i)),'d0'];
    end

    % run the script
    substituteInputParameters(paramToReplace,InputListName,newvalue)    
    
        p = parametersWatercolumn( setupNUMmodel );
        setHTL(mortHTL, mHTL, bHTLquadratic, bHTLdecline)
        setSinkingPOM(p, newSinkingPOM)
        p.tEnd = runningTime;
        sim = simulateWatercolumn(p,lat_to_find,lon_to_find);
        simWCSeason = calcFunctionsTh(sim);
        save(['sim_remin2_',num2str(remin2),'.mat'],'sim');    
        save(['simWC_remin2_',num2str(remin2),'.mat'],'simWCSeason');    
        Ntotal(iparam,:) = simWCSeason.Ntot;
        NPP(iparam,:) = simWCSeason.ProdNet; % (mgC / m2 / day)
        % NPP_annual(i) = simWCSeason.ProdNetAnnual;
end

%%
final_years=3;
monthly_NPP=zeros(length(newRemin2),noYears,12);
monthly_NPP_mean=zeros(length(newRemin2),12);

for iparam = 1:length(newRemin2)
    for i=noYears-3:noYears
         monthly_NPP(iparam,i,:)=reshapeCellToArrayAvg(NPP(iparam,:),i);
    end
  % monthly NPP averaged over the last 3 years of the simulation
    monthly_NPP_mean(iparam,:)=squeeze(mean(monthly_NPP(iparam,noYears-3:noYears,:),2));
end
save(['NPP_monthly_mean_remin2',num2str(remin2),'.mat'],'monthly_NPP_mean');    

%% -------------------------------------------
%                  PLOTS
%-------------------------------------------

figure(1)
clf(1)
tiledlayout(1,length(newRemin2))

nexttile
for i = 1:length(newRemin2)
        plot(monthly_NPP_mean(i,:), 'o-r', 'LineWidth', 1)
        hold on
        plot(NPP_extracted(1,:), 'ob--')
        plot(NPP_extracted(2,:), 'og--')
        plot(NPP_extracted(3,:), 'om--')

legend('NUM_{last3year.avg}','Eppley model', 'Standard VGPM', 'CAFE','Location','best')
xlabel('Time (month)')
ylabel('NPP (mgC / m^2 /day)')
my_title = append('remin2 = ', string(newRemin2(i)));
mTitle = append('Lat: ', string(lat_to_find), ', Lon: ', string(lon_to_find),', ', my_title);
title(mTitle)
end

nexttile
geoscatter(lat_to_find, lon_to_find, 'filled')
exportgraphics(gcf,[saved_png])       
  
saved_png2 = append(name_png,'_details.png');
                
figure(2)
clf(2)
set(gcf,'Position',get(0,'ScreenSize'))
t=tiledlayout(length(newRemin2)/2,length(newRemin2)/2)
for iparam = 1:length(newRemin2)

    nexttile                                   
        plot(squeeze(monthly_NPP_mean(iparam,:)), 'o-r', 'LineWidth', 1)% mgC/m2/day
        hold on
        plot(NPP_extracted(1,:), 'ob--')
        plot(NPP_extracted(2,:), 'og--')
        plot(NPP_extracted(3,:), 'om--')
    xlabel('Time (month)')
    ylabel('NPP (mgC / m^2 /day)')
    mTitle = append('Lat: ', string(lat_to_find), ', Lon: ', string(lon_to_find));
    my_title = append('remin2 = ', string(newRemin2(iparam)));
    title(my_title)
end
title(t,mTitle, 'FontSize', 24)
lg  = legend('NUM','Eppley model', 'Standard VGPM', 'CAFE','NumColumns',2); 
lg.Layout.Tile = 'South'; % <-- Legend placement with tiled layout

exportgraphics(gcf,[saved_png2])       

end
