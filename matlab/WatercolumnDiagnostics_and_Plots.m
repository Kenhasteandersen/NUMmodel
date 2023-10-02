% Water column diagnostic plots
% for graphs it calls function: 
% stackedBarplot_percentage(BpnmVec,legend_names,xlabel_names,mTitle)
%
 lat=0;
 lon=-172;
 p=setupNUMmodel;
 p=parametersWatercolumn(p);
 p.tEnd= 5*365;
 sim=simulateWatercolumn(p,lat,lon);

 % Calculate all photrophic biomass
 Bphyto = calcPhytoplanktonBiomass(sim);

% Calculate Biomass of organism that rely more than 60% on photosynthesis
 simR=getSimRates(sim);
 Bph=calcPhyto(sim,simR);

%% Select day/time and depth layer for plotting
day=p.tEnd-170;
depth_layer=10;
Bphyto4plots = squeeze(mean(Bphyto(:,depth_layer,:),1));
Bph4plots = squeeze(mean(Bph(:,depth_layer,:),1));

Bphyto_pnm=calcPhytoPicoNanoMicroRadius(Bphyto4plots,p,sim);
Bph_pnm=calcPhytoPicoNanoMicroRadius(Bph4plots,p,sim);

nameGroups = sim.p.nameGroup;
typeGroups = sim.p.typeGroups;
Bnum=squeeze(mean(sim.B(:,depth_layer,:),1));
BpnmNUM = calcPicoNanoMicroRadius(Bnum,sim.p,true);
% BpnmNUMinit=BpnmNUM;
BpnmNUMlabels=[BpnmNUM sim.p.typeGroups'] ;

zero_rows= find(all(BpnmNUMlabels(:,1:end-1)==0,2)); % find all-zero rows
nameGroups(zero_rows)=[];
typeGroups(zero_rows)=[];
BpnmNUM(zero_rows, :)=[]; 

% --------------------------------------------------------
% Calculate Unicellular, Multicellular POM Biomass
%--------------------------------------------------------
% Check if there are unicellular groups
if ( any(typeGroups<10) )
   BpnmU=zeros(1,3);
    for iGroup = 1:length(typeGroups)
        if (typeGroups(iGroup)<10)
            BpnmU= BpnmU + BpnmNUM(iGroup,:);
        end
    end
end
iTiles=1;

% Check if there are copepod groups
if ( any(typeGroups>=10) && any(typeGroups<100) )
    iTiles=2;
    BpnmM=zeros(1,3);
     for iGroup = 1:length(typeGroups)
       if (typeGroups(iGroup)==10 || typeGroups(iGroup)==11)
        BpnmM= BpnmM + BpnmNUM(iGroup,:);
       end
    end
end

% Check if there are POM groups
if ( any(typeGroups==100)  )
    BpnmPOM=zeros(1,3);
     for iGroup = 1:length(typeGroups)
       if (typeGroups(iGroup)==100)
        BpnmPOM= BpnmPOM + BpnmNUM(iGroup,:);
       end
    end
end

if (exist('BpnmU','var') && exist('BpnmM','var') && exist('BpnmPOM','var'))
        general_groups = {'Unicellular','Multicellular','POM'};
        ZZ = [  BpnmU; BpnmM; BpnmPOM];

elseif (exist('BpnmU','var') && exist('BpnmM','var'))
        general_groups = {'Unicellular','Multicellular'};
        ZZ = [  BpnmU; BpnmM];

elseif (exist('BpnmM','var') && exist('BpnmPOM','var'))
        general_groups = {'Multicellular','POM'};
        ZZ = [BpnmM; BpnmPOM];
end

%%----------------------------------------------------------------------
% -----------------------------------------------------------------------
%             P L O T S
%
%........................................................................
newcolors = {'#b25781', '#66b9bb','#b5fded',' #b5eded','#b5dded',...
    '#b5cded','#b5bded', '#d7bde2'};
%........................................................................


 for iGroup = 1:length(typeGroups)
    legend_names{iGroup} = [ (nameGroups{iGroup})];
 end

BpnmVec = [Bph_pnm;BpnmU-Bph_pnm; BpnmM;BpnmPOM];
BpnmVec_all = [Bph_pnm;BpnmU ;BpnmU-Bph_pnm; BpnmM; BpnmPOM;BpnmU+BpnmM+BpnmPOM];

tiles = tiledlayout(2,2);
    legend_namesVec=[{'Phyto_{0.6}'},{'U-Phyto_{0.6}'},{'M'},{'POM'}];
    legend_names=[{'Phyto_{0.6}'},{'U'},{'U-Phyto_{0.6}'},{'M'},{'POM'},{'\SigmaB'}];
    size_names={'Pico','Nano','Micro'};

% TOP TILE
nexttile
    stackedBarplot_percentage(Bph_pnm',legend_names,size_names,'Biomass of phytoplankton 60%')
    axis tight

nexttile
    stackedBarplot_percentage(BpnmVec_all,size_names,legend_names,'Biomass of plankton')
    axis tight

nexttile
    stackedBarplot_percentage(BpnmVec',legend_namesVec,size_names,'Biomass of plankton')
    axis tight


nexttile
    bar(ZZ'./sum(ZZ)', 'stacked') 
    legend(general_groups,'Location','bestoutside')
    xticklabels(size_names)
    axis tight
    title(append('Relative biomass of general groups (%) '))

mTitle = append('Biomass partitioning (day: ',num2str(day) ,', depth: ', num2str(sim.z(depth_layer)) ,'m, lat: ',num2str(lat),', lon: ',num2str(lon),')');

title(tiles,mTitle, 'FontSize', 20)

%%
exportgraphics(gcf,[append('biomassNUM','_barplots.png')])       
