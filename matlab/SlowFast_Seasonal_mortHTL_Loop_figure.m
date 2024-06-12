%% Figur for SlowFast article:
% Figure with optimal vulnerability as a function of mortHTL and time

%% Process the files from each mortHTL run
resultfolder=fullfile('ModelResults','seasonal');
thefilestoload=dir(fullfile(resultfolder,'workspace_mortHTL*'));
thefilestoload={thefilestoload.name};
thismortHTL=nan(1,length(thefilestoload));
for fnr=1:length(thefilestoload)
    thismortHTL(fnr)=str2double(extractBetween(thefilestoload(fnr),'mortHTL_0_','.mat'));
    load(fullfile(resultfolder,char(thefilestoload(fnr))));
    Phag_fnr=nan(size(B_phag,1),length(sim.p.nameGroup)-1);
    for i=2:length(sim.p.nameGroup)
    Phag_fnr(:,i-1)=sum(B_phag(:,p.ixStart(i)-2:p.ixEnd(i)-2),2);
    end
[i,j]=max(Phag_fnr');
theoptimal(fnr).group=j;
theoptimal(fnr).time=sim.t(ix);
clear Phag_fnr
end
%% Interpret all model results into same timeline
t=linspace(min(theoptimal(2).time),max(theoptimal(2).time),5000);
optimalgroup=ones(max(thismortHTL),length(t));
for fnr=1:length(thefilestoload)
optimalgroup(thismortHTL(fnr),:)=interp1(theoptimal(fnr).time,theoptimal(fnr).group,t);
end
%% Plot result
figure('Color','w');
time=t(2:end)./365;
mortvector=thismortHTL./10;
palatability_optimal=palatability_these(round(optimalgroup(:,2:end)));
% imagesc(time,mortvector,palatability_optimal);
imagesc(time,mortvector,optimalgroup(:,2:end));
colorbar
xlabel('years since initiation')
ylabel('mHTL')
set(gca,'YDir', 'normal','XGrid','on','GridAlpha',1,'GridLineWidth',1,'XColor','r')

ax=gca;
xt=ax.XTick;
xtminor=min(xt)-1:0.25:max(xt);
ax.XAxis.MinorTickValues=xtminor;
set(gca,'XMinorTick','on','XMinorGrid','on','MinorGridLineStyle','--','MinorGridAlpha',1,'MinorGridLineWidth',0.5)
