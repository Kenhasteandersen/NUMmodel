% Dette script er brugt til at køre NUM med seasonal for forskellige mortHTL.
% Det kører som et loop og laver output "workspace_mortHTL_0_xx.mat hvor xx
% er mortHTL vectoren. Dette output plottes efterfølgende med funktionerne:
% SlowFast_figures_Seasonal_loadandplot.m : surf plot
% number of generalists
n = 20;
% number of generlist groups with different vulnerabilities:
k = 5;
% how many years to run?
tEnd=30;


% Which mortHTL to loop:
morthtlvector=[0.4 0.6 0.7 0.8 0.9];

%place to calculate
latlon=[60,-15];

%year to save
year_plot_start=10;

% Set parameters
[mMinGeneralist,mMaxGeneralist,epsilonL,alphaL,rLstar,alphaN,rNstar,...
    epsilonF,cF,beta,sigma,cLeakage,delta,alphaJ,cR,remin2,reminF]=parametersFastSlow;


%% Calculate vulnerability
% phagotrophic affinity vector
% Get pairs of clearance rate-predation risk: 
alphamin=-4;
alphamax=0;
alphanrintervals=20;
[v,~,alpha,~]=getVulnformulation(alphamin,alphamax,alphanrintervals);
% pick k numbers to use on the clearance rate-predation risk combination
inr=round(linspace(1,length(alpha),k));
alphaF_these=[alpha(inr)];
palatability_these=[v(inr)];

% if specific ones use instead:
% alphaF_these=[alphaF(5) alphaF(15)];
% palatability_these=[palatability(5) palatability(15)];


theOriginalfile='../input/input_generalists_simpleX_ManyGroups.h';
inputFileName='../input/input_generalists_simpleX.h';
[SUCCESS,MESSAGE,MESSAGEID] = copyfile(theOriginalfile, inputFileName);
if ~SUCCESS
    error('error: copying the original file')
end


paramToReplace={'alphaF';'palatability'};
InputListName={'input_generalists_simple';'input_generalists_simple'};
% change to cell array
alphaF_vector=[];
palatability_vector=[];
for i=1:length(alphaF_these)
    % alphaF_this=[num2str(alphaF(i)),'d0 '];
    alphaF_vector=[alphaF_vector num2str(alphaF_these(i)),'d0 '];
    palatability_vector=[palatability_vector num2str(palatability_these(i)),'d0 '];
end
newvalue{1,1}=alphaF_vector;
newvalue{2,1}=palatability_vector;

substituteInputParameters(paramToReplace,InputListName,newvalue,inputFileName)


%update the other parameters with the right amount of groups:
substituteInputParameters('mMinGeneralist','input_generalists_simple',[num2str(k),'*',mMinGeneralist],inputFileName)
substituteInputParameters('mMaxGeneralist','input_generalists_simple',[num2str(k),'*',mMaxGeneralist],inputFileName)

substituteInputParameters('epsilonL','input_generalists_simple',[num2str(k),'*',epsilonL],inputFileName)
substituteInputParameters('alphaL','input_generalists_simple',[num2str(k),'*',alphaL],inputFileName)
substituteInputParameters('rLstar','input_generalists_simple',[num2str(k),'*',rLstar],inputFileName)
substituteInputParameters('alphaN','input_generalists_simple',[num2str(k),'*',alphaN],inputFileName)
substituteInputParameters('rNstar','input_generalists_simple',[num2str(k),'*',rNstar],inputFileName)

substituteInputParameters('epsilonF','input_generalists_simple',[num2str(k),'*',epsilonF],inputFileName)
substituteInputParameters('cF','input_generalists_simple',[num2str(k),'*',cF],inputFileName)
substituteInputParameters('beta','input_generalists_simple',[num2str(k),'*',beta],inputFileName)
substituteInputParameters('sigma','input_generalists_simple',[num2str(k),'*',sigma],inputFileName)

substituteInputParameters('cLeakage','input_generalists_simple',[num2str(k),'*',cLeakage],inputFileName)
substituteInputParameters('delta','input_generalists_simple',[num2str(k),'*',delta],inputFileName)
substituteInputParameters('alphaJ','input_generalists_simple',[num2str(k),'*',alphaJ],inputFileName)
substituteInputParameters('cR','input_generalists_simple',[num2str(k),'*',cR],inputFileName)

substituteInputParameters('remin2','input_generalists_simple',[num2str(k),'*',remin2],inputFileName)
substituteInputParameters('reminF','input_generalists_simple',[num2str(k),'*',reminF],inputFileName)

for jmortHTL=2 %1:length(morthtlvector)
    mortHTL = morthtlvector(jmortHTL);
    mHTL = 1/500^1.5;
    bQuadraticHTL = false;
    bDecliningHTL = false;
    %% Setup NUM
    p = setupGeneralistsSimpleK(n,k, false);
    % p = parametersChemostat(p, 'lat_lon', latlon);
    p = parametersWatercolumn(p);
    p.tEnd = 365.*tEnd;
    setHTL(mortHTL, mHTL, bQuadraticHTL,bDecliningHTL);

    %% Simulatie NUM
    sim = simulateWatercolumn(p, latlon(1),latlon(2));
    % sim = simulateChemostat(p, 'bUnicellularloss', false);

    %% Find the timesteps to use forward 
    ix=find(sim.t>=365*(year_plot_start-1));

    %% get rates and biomass
    jDOC_tmp=zeros(length(ix),size(sim.rates,2),p.n-2);
    jLreal_tmp=jDOC_tmp;
    jFreal_tmp=jDOC_tmp;
    for j=1:length(ix)
        for i=1:size(sim.rates,2)
        jDOC_tmp(j,i,:)=sim.rates(ix(j),i).rates.jDOC';
        jLreal_tmp(j,i,:)=sim.rates(ix(j),i).rates.jLreal';
        jFreal_tmp(j,i,:)=sim.rates(ix(j),i).rates.jFreal';
        end
    end

    B_phag=sim.B(ix,:,:).*(jFreal_tmp./(jLreal_tmp+jDOC_tmp+jFreal_tmp));
    B_photo=sim.B(ix,:,:).*(jLreal_tmp./(jLreal_tmp+jDOC_tmp+jFreal_tmp));
    B_osmo=sim.B(ix,:,:).*(jDOC_tmp./(jLreal_tmp+jDOC_tmp+jFreal_tmp));
    t=sim.t;
   
    save(fullfile('ModelResults','watercolumn','seasonal',['long_workspace_latlon_',num2str(latlon(1)),'_',num2str(latlon(2)),'_mortHTL_0_',num2str(mortHTL*10)]),'-v7.3')
end


