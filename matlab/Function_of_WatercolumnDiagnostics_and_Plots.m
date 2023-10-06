% Water column diagnostic plots
% for graphs it calls function: 
% stackedBarplot_percentage(BpnmVec,legend_names,xlabel_names,mTitle)
%
%  lat=0;
%  lon=-172;
%  p=setupNUMmodel;
%  p=parametersWatercolumn(p);
%  p.tEnd= 5*365;
%  sim=simulateWatercolumn(p,lat,lon);
% day=p.tEnd-170;
% depth_layer=10;


function Function_of_WatercolumnDiagnostics_and_Plots(lat,lon,depth_layer,day,p,sim,simR)

 % Calculate all photrophic biomass
 Bphyto = calcPhytoplanktonBiomass(sim);

% Calculate Biomass of organism that rely more than 60% on photosynthesis
 Bph=calcPhyto(sim,simR);

%% Select day/time and depth layer for plotting
if day==0
   if (depth_layer~=0)
    Bphyto4plots = squeeze(mean(Bphyto(:,depth_layer,:),1));
    Bph4plots = squeeze(mean(Bph(:,depth_layer,:),1));
    Bnum=squeeze(mean(sim.B(:,depth_layer,:),1));
    N=squeeze(mean(sim.N(:,depth_layer),1));
    DOC=squeeze(mean(sim.DOC(:,depth_layer),1));
    Si=squeeze(mean(sim.Si(:,depth_layer),1)); % take care of this if no diatoms

   else
    Bphyto4plots = squeeze( mean(mean(Bphyto(:,:,:),2),1));
    Bph4plots = squeeze(mean(mean(Bph(:,:,:),2),1));
    Bnum=squeeze(mean(mean(sim.B(:,:,:),2),1));
    N=squeeze(mean(mean(sim.N(:,:),2),1));
    DOC=squeeze(mean(mean(sim.DOC(:,:),2),1));
    Si=squeeze(mean(mean(sim.Si(:,:),2),1));
   end
else
   if (depth_layer~=0)
    Bphyto4plots = squeeze(Bphyto(day,depth_layer,:));
    Bph4plots = squeeze(Bph(day,depth_layer,:));
    Bnum=squeeze(sim.B(day,depth_layer,:));
    N=squeeze(sim.N(day,depth_layer));
    DOC=squeeze(sim.DOC(day,depth_layer));
    Si=squeeze(sim.Si(day,depth_layer));
   else
    Bphyto4plots = squeeze( mean(Bphyto(day,:,:),2));
    Bph4plots = squeeze(mean(Bph(day,:,:),2));
    Bnum=squeeze(mean(sim.B(day,:,:),2));
    N=squeeze(mean(sim.N(day,:),2));
    DOC=squeeze(mean(sim.DOC(day,:),2));    
    Si=squeeze(mean(sim.Si(day,:),2));
   end
end

Bphyto_pnm=calcPhytoPicoNanoMicroRadius(Bphyto4plots,p,sim);
Bph_pnm=calcPhytoPicoNanoMicroRadius(Bph4plots,p,sim);

nameGroups = sim.p.nameGroup;
typeGroups = sim.p.typeGroups;
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

%% ----------------------------------------------------------------------
% -----------------------------------------------------------------------
%             P L O T S
%
%........................................................................
newcolors = {'#b25781', '#66b9bb','#b5fded',' #b5eded','#b5dded',...
    '#b5cded','#b5bded', '#d7bde2'};
%........................................................................


 for iGroup = 1:length(typeGroups)
    legendNames{iGroup} = [ (nameGroups{iGroup})];
 end

BpnmVec = [Bph_pnm;BpnmU-Bph_pnm; BpnmM;BpnmPOM];
BpnmVec_all = [Bph_pnm;BpnmU ;BpnmU-Bph_pnm; BpnmM; BpnmPOM;BpnmU+BpnmM+BpnmPOM];
BpnmNUMgeneral = [BpnmU; BpnmM;BpnmPOM];

RR=BpnmNUM./sum(BpnmNUM);
rr=BpnmNUM;
RRR=sum(BpnmNUMgeneral,2);%./sum(BpnmNUMgeneral,"all");
set(gcf, 'Position', get(0, 'Screensize'))
tiles = tiledlayout(3,3);
    legend_namesVec=[{'Phyto_{0.6}'},{'U-Phyto_{0.6}'},{'M'},{'POM'}];
    legend_names=[{'Phyto_{0.6}'},{'U'},{'U-Phyto_{0.6}'},{'M'},{'POM'},{'\SigmaB'}];
    legend_namesNUM=[{'U'},{'M'},{'POM'}];

    size_names={'Pico','Nano','Micro'};

% TOP TILE
nexttile
    stackedBarplot_percentage(Bph_pnm',legend_names,size_names,'Biomass of phytoplankton 60%')
    axis tight

nexttile
    stackedBarplot_percentage(BpnmVec_all,size_names,legend_names,'Biomass of plankton')
    axis tight
 
nexttile
    bar(RR', 'stacked') 
    legend(legendNames,'Location','bestoutside');
    title('Relative biomass of all groups (%) ')
    xticklabels(size_names)
    axis tight
    colororder(newcolors);

nexttile
    stackedBarplot_percentage(sum(BpnmNUMgeneral,2)',legend_namesNUM,"",'Biomass of general plankton groups')
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

legend_namesNUMnn=[{'N'},{'DOC'},{'Si'}];
nexttile % TEST ADDING NUTRIENTS
    stackedBarplot_percentage([N/5.68; DOC; Si/3.4]',... % convert to C-units
        legend_namesNUMnn,"",'NUtrient concentrations')
    axis tight
    
if day==0
   day_str='average';
else
   day_str=num2str(day);
end

 if (depth_layer==0)
    mTitle = append('Biomass partitioning (day: ',day_str ,', depth-integrated, lat: ',num2str(lat),', lon: ',num2str(lon),')');
 else
   mTitle = append('Biomass partitioning (day: ',day_str,', depth: ', num2str(sim.z(depth_layer)) ,'m, lat: ',num2str(lat),', lon: ',num2str(lon),')');
 end
title(tiles,mTitle, 'FontSize', 20)

%%
% exportgraphics(gcf,[append('biomassNUM','_barplots.png')])       
