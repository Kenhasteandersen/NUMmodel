%
% Save the simulation results in an excel file savedSimulaions.xlsx (used in testNewCommit)
% The excel file is in the folder testNewCommit in the current folder
%
% For chemostat: save N, DOC, Si, total B, jLreal, jFreal & jDOC at the end of the simulation
% For watercolumn: save the average over the 100 first meters of N, DOC, Si,
% total B, jLreal, jFreal & jDOC at the end of the simulation
% For global: save the average over the 100 first meters of N, DOC, Si,
% total B, jLreal, jFreal & jDOC at the end of the simulation on 4 different zones
%
% In:
%   sim - simulation structure
%   rates - jLreal, jFreal & jDOC (use calcRates)
%   comment - a comment on the simulation, empty by delauft (must be string)
%

function saveResultsSimulation(sim,rates,comment)

arguments
    sim struct
    rates struct
    comment string = ''
end

%
% Calculates Data
%
switch sim.p.nameModel
    case 'chemostat'
        N = sim.N(end);
        DOC = sim.DOC(end);
        Si = sim.Si(end);
        B = sum(sim.B(end,:));
        jLreal = sum(rates.jLreal(end,:));
        jFreal = sum(rates.jFreal(end,:));
        jDOC = sum(rates.jDOC(end,:));
        U = {N DOC Si B jLreal jFreal jDOC};

    case 'watercolumn'
        N = mean(sim.N(end,sim.z<=100));
        DOC = mean(sim.DOC(end,sim.z<=100));
        Si = mean(sim.Si(end,sim.z<=100));
        B = mean(sum(sim.B(end,sim.z<=100,:),3));
        jLreal = mean(sum(rates.jLreal(end,sim.z<=100,:),3));
        jFreal = mean(sum(rates.jFreal(end,sim.z<=100,:),3));
        jDOC = mean(sum(rates.jDOC(end,sim.z<=100,:),3));
        U = {N DOC Si B jLreal jFreal jDOC};

    case 'global'
        %
        % Zones
        %
        load('testNewCommit\zones.mat','zones'); %indo_pacific,atlantic,southern,northern
        % plots zones
        bathy = squeeze(sim.bathy(:,:,1));
        figure 
        hold on
        for i=1:length(zones)
            data = bathy.*zones{i}*i; data(zones{i}==0)=NaN;
            contourf(sim.x,sim.y,data',5,'linestyle','none')
        end
        contour(sim.x,sim.y,isnan(squeeze(sim.N(end,:,:,1))'),1,'LineWidth',3.5)
        hold off
        colorbar('Ticks',0:length(zones),...
         'TickLabels',{'Land','Indio-Pacific Warm Water','Atlantic Warm Water','Southern Cold Water','Northen Cold Water'})
        colormap([0 0 0; 1 0.25 0.08; 1 0.95 0.12; 0.4 0.4 0.27; 0.5 0.5 0.5]);
        clim([0 4])
        xlabel('longitude (°)')
        ylabel('latitude (°)')
        title('Zones')

        %datas
        dataN = squeeze(sim.N(end,:,:,sim.z<=100)); 
        dataDOC = squeeze(sim.DOC(end,:,:,sim.z<=100));
        dataSi = squeeze(sim.Si(end,:,:,sim.z<=100));
        dataB = squeeze(sum(sim.B(end,:,:,sim.z<=100,:),5));
        dataJLreal = squeeze(sum(rates.jLreal(:,:,sim.z<=100,:),4));
        dataJFreal = squeeze(sum(rates.jFreal(:,:,sim.z<=100,:),4));
        dataJDOC = squeeze(sum(rates.jDOC(:,:,sim.z<=100,:),4));

        data = {dataN dataDOC dataSi dataB dataJLreal dataJFreal dataJDOC};

        % Mean in surface layer (100m) in each zones
        U = cell(1,length(zones)*length(data));
        for i = 1:length(zones) 
            for j = 1:length(data)
                Data = data{j}(zones{i}==1);
                U{j+length(data)*(i-1)} = mean(Data(~isnan(Data)),'all');
            end
        end
end
U=[{datetime, comment},U];

%
% Saves Data
%
if ~exist('testNewCommit\savedSimulations.xlsx','file') %check if the excel file already exists
    %creates the table for chemostat and watercolumn
    varTypes = ["datetime","string","double","double","double","double","double","double","double"];
    varNames = ["Date/Time","Comment","N","DOC","Si","Biomass","jLreal","jFreal","jDOC"];
    Size = [0 length(varNames)];
    T1 = table('Size',Size,'VariableTypes',varTypes,'VariableNames',varNames);
    T1 = [{datetime,'creation'}, cell(1,Size(2)-2);T1];
    %creates the table for global
    varTypes = ["datetime","string","double","double","double","double","double","double","double","double","double","double","double","double","double","double","double","double","double","double","double","double","double","double","double","double","double","double","double","double"];
    varNames = ["Date/Time","Comment","N_indoPacific","DOC_indoPacific","Si_indoPacific","Biomass_indoPacific","jLreal_indoPacific","jFreal_indoPacific","jDOC_indoPacific","N_atlantic","DOC_atlantic","Si_atlantic","Biomass_atlantic","jLreal_atlantic","jFreal_atlantic","jDOC_atlantic","N_south","DOC_south","Si_south","Biomass_south","jLreal_south","jFreal_south","jDOC_south","N_north","DOC_north","Si_north","Biomass_north","jLreal_north","jFreal_north","jDOC_north"];
    Size = [0 length(varNames)];
    T2 = table('Size',Size,'VariableTypes',varTypes,'VariableNames',varNames);
    T2 = [{datetime,'creation'}, cell(1,Size(2)-2);T2];
    
    %creates the excel file, one sheet for each model
    writetable(T1,'testNewCommit\savedSimulations.xlsx','FileType','spreadsheet','Sheet','chemostat')
    writetable(T1,'testNewCommit\savedSimulations.xlsx','FileType','spreadsheet','Sheet','watercolumn')
    writetable(T2,'testNewCommit\savedSimulations.xlsx','FileType','spreadsheet','Sheet','global')
    
    %saves the new data in the file
    switch sim.p.nameModel
        case 'global'
            T = T2;
        otherwise
            T = T1;
    end
else
    %loads the sheet 
    T=readtable('testNewCommit\savedSimulations.xlsx','FileType','spreadsheet','Sheet',sim.p.nameModel,'VariableNamingRule','preserve');
end

T = [U;T];
writetable(T,'testNewCommit\savedSimulations.xlsx','FileType','spreadsheet','Sheet',sim.p.nameModel)