% Returns stacked barplots of pico-, nano- and microplankton
% for different type groups
% example of Bnum: Bnum=squeeze(mean(simNUM.B(:,1,:),1));

function [tiles] = PicoNanoMicroBarplots(sim,Bnum) 

nameGroups = sim.p.nameGroup;
typeGroups = sim.p.typeGroups;
BpnmNUM = calcPicoNanoMicroRadius(Bnum,sim.p,true);
% BpnmNUMinit=BpnmNUM;
BpnmNUMlabels=[BpnmNUM sim.p.typeGroups'] ;

zero_rows= find(all(BpnmNUMlabels(:,1:end-1)==0,2)); % find all-zero rows
nameGroups(zero_rows)=[];
typeGroups(zero_rows)=[];
BpnmNUM(zero_rows, :)=[]; 

%--------------------------------------------------------
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
%........................................................................
newcolors = {'#b25781', '#66b9bb','#b5fded',' #b5eded','#b5dded',...
    '#b5cded','#b5bded', '#d7bde2'};
%........................................................................

RR=BpnmNUM./sum(BpnmNUM);
rr=BpnmNUM;

 for iGroup = 1:length(typeGroups)
    legend_names{iGroup} = [ (nameGroups{iGroup})];
 end

clf
tiles = tiledlayout(iTiles,2);

nexttile
bar(RR', 'stacked') 
legend(legend_names,'Location','bestoutside');
title('Relative biomass of all groups (%) ')
xticklabels({'Pico','Nano','Micro'})
axis tight
colororder(newcolors);

nexttile
bar(rr', 'stacked') 
ylabel('Biomass (\mugC/L)')
xticklabels({'Pico','Nano','Micro'})
axis tight
title('Absolute biomass of all groups')

if iTiles==2
    nexttile
    bar(ZZ'./sum(ZZ)', 'stacked') 
    legend(general_groups,'Location','bestoutside')
    xticklabels({'Pico','Nano','Micro'})
    axis tight
    title(append('Relative biomass of general groups (%) '))

    nexttile
    bar(ZZ', 'stacked') 
    ylabel('Biomass (\mugC/L)')
    xticklabels({'Pico','Nano','Micro'})
    axis tight
    title('Absolute biomass of general groups')
end

