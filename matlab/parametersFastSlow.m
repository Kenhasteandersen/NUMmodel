function [mMinGeneralist,mMaxGeneralist,epsilonL,alphaL,rLstar,...
    alphaN,rNstar,epsilonF,cF,beta,sigma,cLeakage,delta,alphaJ,cR,remin2,reminF]=parametersFastSlow()
% Denne funktion kalder parametrene der benyttes i dette studie og som
% sikrer at der ikke er brugt forskellige parametre undervejs

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
end
