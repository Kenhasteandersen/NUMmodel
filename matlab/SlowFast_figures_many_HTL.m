mAdult=[];

%% define parameters
% number of size groups:
n = 20;
% number of generlist groups:
k = 10;
% run in parallel?:
bParallel = false;
% Years of simulation
tEnd=5;
% Years to save;
tSave=2;
% nutrients
d=0.1;
% set HTL parameters:
%the HTL mortality (or the constant if bQuadratic=true)
mortHTL=0.1;
mHTL=1/500^1.5;
bQuadraticHTL=false;        
bDecliningHTL=false;

%calculate 
alphaF=logspace(-4,-2,10);
meanval=0.8373;
% vulnera=(10.^(P(2)+P(1).*log10(alphaF)))./meanval;


%% setup simulation
% % % % p = setupGeneralistsSimpleK(n,k, bParallel);
% % % % p = parametersChemostat(p);
% % % % p.tEnd = 365*tEnd;
% % % % p.d = d;
% % % % 
% % % % % Change colors for groups of Generalists:
% % % % r=[13 15 56 119 189];
% % % % g=[139 178 192 207 227];
% % % % b=[217 242 242 242 242];
% % % % j=linspace(0.3,0.9,k);
% % % % 
% % % % blue(1,:)=[58/255,67/255,186/255];
% % % % blue(2,:)=[4/255,146/255,194/255];
% % % % blue(3,:)=[130/255,237/255,253/255];
% % % % for i=4:k+1
% % % %  blue(i,:)=[130/255,237/255,253/255];
% % % % end
% % % % 
% % % % for i=2:k+1
% % % %     %p.colGroup{i}=[0,0,j(i-1)];
% % % %     p.colGroup{i}=blue(i-1,:);
% % % % end
% % % % %Change color of N so I can see generalists
% % % % 
% % % % p.colNutrients{1}=[0,1,1];
% % % % 
% % % % %p.umin=p.umin.*0.5;
% % % % %p.u0(p.idxB:p.n) = 0.5;
% % % % %
% % % % % Set to "normal" HTL mortality if there are no copepods:
% % % % %
% % % % 
% % % % setHTL(mortHTL, mHTL, bQuadraticHTL,bDecliningHTL);

%% Simulate
usesaved='yes';
mortHTL_vector=linspace(0,0.3,20);
dvector=mortHTL_vector;
switch usesaved
    case 'yes'
        load('bio_saved.mat')
    case 'no'
        % iterate over d=0.001 to 0.1
        
        B=nan(length(dvector),n*(k+1));
        jDOC=B;
        jLreal=B;
        jFreal=B;
        
        for i=1:length(mortHTL_vector)
            setHTL(mortHTL_vector(i), mHTL, bQuadraticHTL,bDecliningHTL);
            sim = simulateChemostat(p, 100);
            % create average of last
            tsave_start=365*(tEnd-tSave)+1;
            idx=find(sim.t>tsave_start);

            B(i,:)=mean(sim.B(idx,:),1,'omitnan');
            jDOC_tmp=nan(sim.p.n-2,length(idx));
            jLreal_tmp=jDOC_tmp;
            jFreal_tmp=jDOC_tmp;
            for j=1:length(idx)
                jDOC_tmp(:,j)=sim.rates(idx(j)).jDOC;
                jLreal_tmp(:,j)=sim.rates(idx(j)).jLreal;
                jFreal_tmp(:,j)=sim.rates(idx(j)).jFreal;
            end
            jDOC(i,:)=mean(jDOC_tmp,2,'omitnan')';
            jLreal(i,:)=mean(jLreal_tmp,2,'omitnan')';
            jFreal(i,:)=mean(jFreal_tmp,2,'omitnan')';
        end
end



%% Plot
% standard plots:
load('cmaps.mat')
B_DOC=B.*(jDOC./(jDOC+jLreal+jFreal));
B_L=B.*(jLreal./(jDOC+jLreal+jFreal));
B_F=B.*(jFreal./(jDOC+jLreal+jFreal));
themax=max([max(max(B_DOC)) max(max(B_L)) max(max(B_F))]);
%%
PlotNutrientVSSize(p,B,'jet','total biomass',max(max(B)))
%%
PlotNutrientVSSize(p,B_DOC,cmap_osmotrophy,'osmotrophic biomass',themax)
%%
PlotNutrientVSSize(p,B_L,cmap_phototrophy,'phototrophic biomass',themax)
%%
PlotNutrientVSSize(p,B_F,cmap_phagotrophy,'phagotrophic biomass',themax)
%% plot biomassVSnutrint divided into different groups
figure('color','w');
subplot(3,1,1);
%plot(dvector,sum(B(:,1:20),2),'b-');
hold on; 
plot(dvector,sum(B(:,21:40),2),'g-');
set(gca,'XScale','log');
box on
ylabel('biomass (\mug C l^{(-1)})','FontSize', 14)


subplot(3,1,2);
plot(dvector,sum(B(:,41:50),2),'r-');
set(gca,'XScale','log');
ylabel('biomass (\mug C l^{(-1)})','FontSize', 14)


subplot(3,1,3);
plot(dvector,sum(B(:,51:end),2),'k-')
set(gca,'XScale','log')
xlabel('mixing rate (day^{-1})','FontSize', 14)
ylabel('biomass (\mug C l^{(-1)})','FontSize', 14)

%%
thism=p.m(23:42);
gnr=1:10;
[x2,y2]=meshgrid(thism,gnr);
figure('color','w');
for dd=1:size(B,1)
for i=1:10
    grupnr=i+1;
    thisB(i,:)=B(dd,p.ixStart(grupnr)-2:p.ixEnd(grupnr)-2);
end
nexttile
surf(x2,y2,thisB)
set(gca,'XScale','log','YScale','log');xlabel('cell mass');ylabel('generalist group');zlabel('biomass')
view(gca,[31.5832061068702 22.2]);
subtitle(['mortHTL ',num2str(mortHTL_vector(dd))]);
zlim([0 50])
end



%%
nexttile(5)
B1=B(:,p.ixStart(2)-2:p.ixEnd(2)-2);
%B1(B1==0)=0.000001;
B2=B(:,p.ixStart(3)-2:p.ixEnd(3)-2);
%B2(B2==0)=0.000001;
Bratio=B1./(B1+B2);
% Bratio(Bratio<=1)=NaN;
% Bratio2=B2./B1;
% Bratio2(Bratio2<=1)=NaN;
pcolor(x,y,Bratio);

shading flat
set(gca,'YTick', [])
set(gca,'XScale','log','YScale','log')
theXspace=10.^(linspace(log10(p.m(p.ixStart(i))),log10(p.m(p.ixEnd(i))),4));
ax=gca;ax.XTick=theXspace;
ax=gca;ax.XAxis.TickLabels = compose('%g', ax.XAxis.TickValues);
ax.FontSize = 16;
xlabel('Cell mass (\mug C)','FontSize', 20)

title('Generlist 2/generalist 1')
colorbar



%% plotPHTL
% this plot is made just for my own understanding of the mortality function
% can be turned on if wanted to but it is also in the word document
%figure('color','w');
%plot(m(1:10),1./(1+(m(1:10)./(1/500^1.5)).^(-2)),'-*')
%set(gca,'XScale','log','YScale','linear')
%% Plot functions
function PlotNutrientVSSize(p,Bio,cmap,figuretitle,themax)
figure('color','w','Name',figuretitle);
if p.nGroups==3
    tiledlayout(2,p.nGroups+1,'TileSpacing','compact')
end

theX=logspace(-7,2,10);
%dvector=logspace(-3,-1,20);
dvector=linspace(0.05,0.5,10);
for i=1:p.nGroups
    nexttile(i)
    %subplot(1,p.nGroups,i)
    m=p.m(p.ixStart(i):p.ixEnd(i));
    [x,y]=meshgrid(p.m(p.ixStart(i):p.ixEnd(i)),dvector);
    contourf(x,y,Bio(:,p.ixStart(i)-2:p.ixEnd(i)-2),100,'LineStyle','none');
    shading flat
    hold on
    thexlim=xlim; theylim=ylim;
    plot([thexlim(1) thexlim(1) thexlim(2) thexlim(2) thexlim(1)],[theylim(1) theylim(2) theylim(2) theylim(1) theylim(1)],'-k')
    if i~=1;set(gca,'YTick', []);end
    set(gca,'XScale','log','FontSize',14)
    ax=gca;
    test=theX(theX(1,:) >= ax.XLim(1) & theX(1,:) <= ax.XLim(2));
    %ax.XTick=[0.0000001 0.00001 0.001 0.1];
    ax.XTick=test;
    % theXspace=10.^(linspace(log10(p.m(p.ixStart(i))),log10(p.m(p.ixEnd(i))),4));
    % ax=gca;ax.XTick=theXspace;
    % ax=gca;ax.XAxis.TickLabels = compose('%g', ax.XAxis.TickValues);
    xlabel('Cell mass (\mug C)','FontSize', 14)
    colormap(cmap);    
    if i==1;ylabel('mortHTL','FontSize', 14);end
    title([p.nameGroup{i},' ',num2str(length(find(p.typeGroups(1:i)==p.typeGroups(i))))],'FontSize', 16,'FontWeight','normal');
    clim([0 themax])
    box on
    %if i==3;colorbar;end
     
end
nexttile(p.nGroups+1)
colormap(cmap)
cbar=colorbar('Location','westoutside');
ylabel(cbar, '\mug C l^{(-1)}','FontSize', 14)
clim([0 themax])
set(gca,'Visible','off')
end


