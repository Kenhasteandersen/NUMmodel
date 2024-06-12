% This script calculates the "optimal" vulnerability for the phagotrophs in
% a simulation with k numbers of groups.

%% Define Parameters
% number of size groups:
n = 10;
% number of generlist groups with different vulnerabilities:
k = 1;

% HTL mortality
% - HTL mortality pressure (vector)
mortHTLvector=linspace(0.1,1,20);
mortHTL=mortHTLvector(end);
% - HTL lower limit
mHTL=1/500^1.5;
% - HTL use quadratic mortality?
bQuadraticHTL=false;
bDecliningHTL=false;

% mixing rate vector:
dvector=logspace(-3,-1,20);
d=dvector(end);

% phagotrophic affinity vector
alphaF=logspace(-4,0,k);
% parameters that define the vulnerability
gamma=0.6505;
v0_coefficient=1.7198;

% Years of simulation
tEnd=10;
% Make a mean of the last tSave years:
tSave=5;

% run simulation in parallel?:
bParallel = false;
% use already saved simulation?
usesaved='no';

mMinGeneralist = '2.5d-8';  % Smallest cell size [mug C]
mMaxGeneralist = '10.0d0';  %Largest cell size [mug C]

% Light uptake:
epsilonL = '0.8d0';         % Light uptake efficiency []
alphaL = '0.3d0';           % Light affinity coef. [1/(uE/m2/s) 1/day um]
rLstar = '7.5d0';           % Light affinity cross-over size [um]

% Dissolved nutrient and DOC uptake:
alphaN = '0.972d0';         % Diffusive affinity coefficient [L/d/mugC/um^2]
rNstar = '0.4d0';           % Diffusive affinity cross-over size [um]

% Phagotrophy:

epsilonF = '0.8d0';	        % Food assimilation efficiency [-]
cF = '30.0d0';              % Max phagotrophy coefficient [um/day]
beta = '500.d0';            % Preferred predator-prey mass ratio
sigma = '1.3d0';            % Preferred predator-prey mass range

% Metabolism:

cLeakage = '0.03d0';        % Passive leakage of C and N
delta = '0.05d0';           % Thickness of cell wall [um]
alphaJ = '1.5d0';           % Constant for jMax [day-1]
cR = '0.1d0';               % Basal metabolism relative to jMax [-]

% Biogeo:
remin2 = '0.5d0';           % Fraction of viral lysis remineralized to DOC
reminF = '0.1d0';           % Fraction of feeding losses remineralized

%% Calculate vulnerability
palatability=10^(v0_coefficient).*alphaF.^(gamma);
palatability=palatability./mean(palatability);

switch usesaved
    case 'yes'
        load('bio_saved.mat')
    case 'no'

        %% Make input file ready
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
        for i=1:length(alphaF)
            % alphaF_this=[num2str(alphaF(i)),'d0 '];

            alphaF_vector=[alphaF_vector num2str(alphaF(i)),'d0 '];
            palatability_vector=[palatability_vector num2str(palatability(i)),'d0 '];
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


        %% setup simulation
        p = setupGeneralistsSimpleK(n,k, bParallel);
        
        p = parametersChemostat(p,'seasonalAmplitude', 10); % set seasonal options
        % p = parametersChemostat(p);
        % p = parametersChemostat(p,'lat_lon', [60 -10]); % set seasonal options
        p.tEnd = 365*tEnd;
        p.tSave2 = tSave;

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

        %% Simulate
        B=nan(length(d),length(mortHTL),n*(k+1),length(mortHTL));
        jDOC=B;
        jLreal=B;
        jFreal=B;
        % iterate over HTL mortality
        for imort=1:length(mortHTL)
            disp(['loop imort =',num2str(imort),' out of ',num2str(length(mortHTLvector))])
            setHTL(mortHTL, mHTL, bQuadraticHTL,bDecliningHTL);
            % setHTL(mortHTLvector(imort), mHTL, bQuadraticHTL,bDecliningHTL);
            % iterate over mixing rate
            for i=1:length(d)
                disp(['   loop dvector =',num2str(i),' out of ',num2str(length(dvector))])
                tic
                % p.d=d;
                % p.d=dvector(i);
                sim = simulateChemostat(p, 100);
                % create average of result since tSave
                tsave_start=365*(tEnd-tSave)+1;
                idx=find(sim.t>tsave_start);

                B(i,imort,:,:)=sim.B(idx,:);
                jDOC_tmp=nan(sim.p.n-2,length(idx));
                jLreal_tmp=jDOC_tmp;
                jFreal_tmp=jDOC_tmp;
                for j=1:length(idx)
                    jDOC_tmp(:,j)=sim.rates(idx(j)).jDOC;
                    jLreal_tmp(:,j)=sim.rates(idx(j)).jLreal;
                    jFreal_tmp(:,j)=sim.rates(idx(j)).jFreal;
                end
                jDOC(i,imort,:,:)=mean(jDOC_tmp,2,'omitnan')';
                jLreal(i,imort,:,:)=mean(jLreal_tmp,2,'omitnan')';
                jFreal(i,imort,:,:)=mean(jFreal_tmp,2,'omitnan')';
                toc
            end
        end
end



%% Plot
% standard plots:
load('cmaps.mat')
B_DOC=B.*(jDOC./(jDOC+jLreal+jFreal));
B_L=B.*(jLreal./(jDOC+jLreal+jFreal));
B_F=B.*(jFreal./(jDOC+jLreal+jFreal));

%%
% go through each of the nutrient levels (the y direction) and make a
% B_F(vuln.level,size group):
for mortgroup=1:length(mortHTLvector)
    for ngroup=1:length(dvector)
        a=nan(k,n);
        for i=2:k+1
            a(i-1,:)=B_F(ngroup,sim.p.ixStart(i)-2:sim.p.ixEnd(i)-2,mortgroup);
        end


        figure;
        contourf(a)
        set(gca,'XScale','log')
        ylabel('vulnerability group number')
        xlabel('size')

        [ix_max_vuln_group,ix_max_size]=find(a == max(max(a)));
        best_vuln(mortgroup,ngroup)=ix_max_vuln_group;
        best_size(mortgroup,ngroup)=ix_max_size;
    end
end








