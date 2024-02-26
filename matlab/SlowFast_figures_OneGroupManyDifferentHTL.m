mAdult=[];

%% define parameters
% number of size groups:
n = 20;
% number of generlist groups:
k = 1;
% run in parallel?:
bParallel = false;
% Years of simulation
tEnd=5;
% Years to save;
tSave=2;
% timestep:
d=0.1;
% set HTL parameters:
%the HTL mortality (or the constant if bQuadratic=true)
mortHTL=0.1;
mHTL=1/500^1.5;
bQuadraticHTL=false;
bDecliningHTL=false;


% v0_coefficient og gamma fås fra nielsen_kiørboe_2021.m;
gamma=0.6505;
v0_coefficient=1.7198;
%x-aksen på thomas plot
alpha=logspace(-4,0,1000);

%calculate 
% alphaF=logspace(-4,-2,10);
% meanval=0.8373;
% vulnera=(10.^(P(2)+P(1).*log10(alphaF)))./meanval;

% Load the parameters
param_vector=load('alphaF_palatability.mat');

paramToReplace={'alphaF';'palatability'};
% which input list do they belong to?
InputListName={'input_generalists_simple';'input_generalists_simple'};
inputFileName='../input/input_generalists_simpleX_1group.h';
copyFileName='../input/input_generalists_simpleX.h';
theOriginalfile='../input/input_generalists_simpleX_1group_Startup.h';





%% Modify input file
[SUCCESS,MESSAGE,MESSAGEID] = copyfile(theOriginalfile, inputFileName);
if ~SUCCESS
    error('error: copying the original file')
end




%Calculate the vulnerability:
v=10^(v0_coefficient).*alpha.^(gamma);


%% setup simulation
%p = setitup(n,k,bParallel,tEnd,d,mortHTL, mHTL, bQuadraticHTL,bDecliningHTL)


%% Simulate
usesaved='no';
yvector=param_vector.alphaF;
%epsilonL_vector=0.1:0.1:15;
%epsilonL_vector=linspace(0.4,15,10);
switch usesaved
    case 'yes'
        load('bio_saved.mat')
    case 'no'
        % iterate over d=0.001 to 0.1
        
        B=nan(length(yvector),n*(k+1));
        jDOC=B;
        jLreal=B;
        jFreal=B;
        for i=1:length(yvector)
            alphaF=param_vector.alphaF(i);
            % alphaF=epsilonL_vector(i);
            palatability=param_vector.palatability(i);
            %%
            ChangeTHeParam(paramToReplace,InputListName,inputFileName,alphaF,palatability)
            %%
            [SUCCESS,MESSAGE,MESSAGEID] = copyfile(inputFileName, copyFileName);
            p = setitup(n,k,bParallel,tEnd,d,mortHTL, mHTL, bQuadraticHTL,bDecliningHTL);
            p.d=d;
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
figure('color','w');
tiledlayout(4,1,'TileSpacing','compact')

PlotNutrientVSSize(p,B,'jet','total biomass',max(max(B)),yvector,i1)
%%
PlotNutrientVSSize(p,B_DOC,cmap_osmotrophy,'osmotrophic biomass',themax,yvector,i2)
%%
PlotNutrientVSSize(p,B_L,cmap_phototrophy,'phototrophic biomass',themax,yvector,i3)
%%
PlotNutrientVSSize(p,B_F,cmap_phagotrophy,'phagotrophic biomass',themax,yvector,i4)
%% plot biomassVSnutrint divided into different groups
figure('color','w');

subplot(3,1,1);
%plot(dvector,sum(B(:,1:20),2),'b-');
hold on; 
plot(yvector,sum(B(:,21:40),2),'g-');
set(gca,'XScale','log');
box on
ylabel('biomass (\mug C l^{(-1)})','FontSize', 14)


subplot(3,1,2);
plot(yvector,sum(B(:,41:50),2),'r-');
set(gca,'XScale','log');
ylabel('biomass (\mug C l^{(-1)})','FontSize', 14)


subplot(3,1,3);
plot(yvector,sum(B(:,51:end),2),'k-')
set(gca,'XScale','log')
xlabel('mixing rate (day^{-1})','FontSize', 14)
ylabel('biomass (\mug C l^{(-1)})','FontSize', 14)


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
function PlotNutrientVSSize(p,Bio,cmap,figuretitle,themax,dvector,inr)
% figure('color','w','Name',figuretitle);
% if p.nGroups==3
%     tiledlayout(2,p.nGroups+1,'TileSpacing','compact')
% end

theX=logspace(-7,2,10);

for i=2:p.nGroups
    it=nexttile();
    %subplot(1,p.nGroups,i)
    m=p.m(p.ixStart(i):p.ixEnd(i));
    [x,y]=meshgrid(p.m(p.ixStart(i):p.ixEnd(i)),dvector);
    contourf(x,y,Bio(:,p.ixStart(i)-2:p.ixEnd(i)-2),100,'LineStyle','none');
    shading flat
    hold on
    thexlim=xlim; theylim=ylim;
    plot([thexlim(1) thexlim(1) thexlim(2) thexlim(2) thexlim(1)],[theylim(1) theylim(2) theylim(2) theylim(1) theylim(1)],'-k')
    if i~=1;set(gca,'YTick', []);end
    set(gca,'XScale','log','YScale','log','FontSize',14)
    ax=gca;
    test=theX(theX(1,:) >= ax.XLim(1) & theX(1,:) <= ax.XLim(2));
    %ax.XTick=[0.0000001 0.00001 0.001 0.1];
    ax.XTick=test;
    % theXspace=10.^(linspace(log10(p.m(p.ixStart(i))),log10(p.m(p.ixEnd(i))),4));
    % ax=gca;ax.XTick=theXspace;
    % ax=gca;ax.XAxis.TickLabels = compose('%g', ax.XAxis.TickValues);
    xlabel('Cell mass (\mug C)','FontSize', 14)
    colormap(cmap);    
    if i==2;ylabel('strategy type','FontSize', 14);end
    title([p.nameGroup{i},' ',num2str(length(find(p.typeGroups(1:i)==p.typeGroups(i))))],'FontSize', 16,'FontWeight','normal');
    clim([0 themax])
    box on
    %if i==3;colorbar;end
     
end
% nexttile(p.nGroups+1)
colormap(cmap)
cbar=colorbar('Location','westoutside');
ylabel(cbar, '\mug C l^{(-1)}','FontSize', 14)
clim([0 themax])
set(gca,'Visible','off')
end

function ChangeTHeParam(paramToReplace,InputListName,inputFileName,alphaF,palatability)

% change to cell array
combinedvalues=[alphaF,palatability];
newvalue=cell(length(combinedvalues),1);
for i=1:length(combinedvalues)
    newvalue{i}=[num2str(combinedvalues(i)),'d0'];
end

% run the script
substituteInputParameters(paramToReplace,InputListName,newvalue,inputFileName)
end

function p = setitup(n,k,bParallel,tEnd,d,mortHTL, mHTL, bQuadraticHTL,bDecliningHTL)
p = setupGeneralistsSimpleK(n,k, bParallel);
p = parametersChemostat(p);
p.tEnd = 365*tEnd;
p.d = d;

% Change colors for groups of Generalists:
r=[13 15 56 119 189];
g=[139 178 192 207 227];
b=[217 242 242 242 242];
j=linspace(0.3,0.9,k);

blue(1,:)=[58/255,67/255,186/255];
blue(2,:)=[4/255,146/255,194/255];
blue(3,:)=[130/255,237/255,253/255];
for i=4:k+1
 blue(i,:)=[130/255,237/255,253/255];
end

for i=2:k+1
    %p.colGroup{i}=[0,0,j(i-1)];
    p.colGroup{i}=blue(i-1,:);
end
%Change color of N so I can see generalists

p.colNutrients{1}=[0,1,1];

%p.umin=p.umin.*0.5;
%p.u0(p.idxB:p.n) = 0.5;
%
% Set to "normal" HTL mortality if there are no copepods:
%

setHTL(mortHTL, mHTL, bQuadraticHTL,bDecliningHTL);
end