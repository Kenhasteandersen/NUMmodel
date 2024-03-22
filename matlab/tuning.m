function [y,fval,exitflag,output] = tuning(parinit, siminit)
arguments
    parinit = [0.067, 14.5, 0.007]; % kw, u, and mortHTL
    siminit = [];
end

%p = setupNUMmodel(bParallel= true);
p = setupNUMmodel([0.5 2], [10 1000], 6,6,1, bParallel= true); % A fast version of the NUM setup
p = parametersGlobal(p,1);
p.BC_POMclosed = false;
p.BCdiffusion = [1 0 10]; % Diffusivities of N, DOC, and Si at the bottom
setSinkingPOM(p, parinit(2))
p.tEnd = 3*365;
p.bUse_parday_light = true; % Use the PARday file for light
p.kw = parinit(1);

mHTL = .1;
mortHTL = parinit(3);
bHTLdecline = false;
bHTLquadratic = true;

%
% Run a long initial condition:
%
if isempty( siminit )
    setHTL(mortHTL, mHTL, bHTLquadratic, bHTLdecline)
    p.tEnd = 10*365;
    %load('NUMmodelsmall.mat','sim');
    siminit = simulateGlobal(p,bCalcAnnualAverages=true);
    save("NUMmodel.mat",'siminit','-v7.3')
end

%p.tSave = 1;
%p.tEnd = 12;

% The objective we are optimizing towards:

objExpected = [1 1 1 800 200 1000]; % Pico, POC, copepods; NPP , eutrophic oligotrophic and seasonal

options = optimoptions('fmincon');
options.Display = 'iter';

% Min and max vlues of the parameters:
parmin = [0.04, 5, 0.001];
parmax = [0.12, 20, 0.1];

% Run optimization
[y,fval,exitflag,output] = fmincon(@objective, parinit, [],[],[],[], parmin,parmax, [], options);

    function err = objective(par)
        % Set parameters:
        p.kw = par(1);
        setSinkingPOM(p,par(2));
        setHTL(par(3), mHTL, bHTLquadratic, bHTLdecline);

        % Simulate:
        sim = simulateGlobal(p, siminit,bCalcAnnualAverages=true, bVerbose=false);

        % Plots:
        figure(1)
        obj = EvaluateRun(sim);
        figure(2)%f
        plotWatercolumnTime(sim,60,-15,depthMax=200)
        drawnow

        % Evaluate objective criterion:
        obj = obj ./ objExpected;
        obj(4:6) = log(obj(4:6));
        err = double(mean(abs(obj)));
        fprintf("Pars: %f, %f, %f; err: %f\n",[par,err]);
        fprintf("Objective vector: %f, %f, %f, %f, %f, %f\n", obj)
        fprintf("====================\n");
    end

end