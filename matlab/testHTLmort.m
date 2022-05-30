function sim = testHTLmort(mortHTL,mAdult,options,mortm)
%tests the response of NUM model in custom HTL values
%useful for creating bifurcation diagrams HTL-biomass
% for instructions of how to implement the custom library see readme in git

%Input: 
% mAdult: specifies the number of zooplankton groups and their maximum size
% _______________________default: 1 group with maximum size of 20μgC
% mortHTL: HTL mortality per day - default: 0.1
% options: under development, specifies the shape of the HTLmort curve
% 1: constant mortality for all size classes | 2: sigmoid curve
% 3: coming soon          | default : 1
% mortm: mass (μg C) where mortality = 1/2 mortHTL

%Output
% sim: structure with NUM results
arguments
    mortHTL double = 0.1
    mAdult double = [20]
    options double = 1
    mortm double = min(mAdult)
end

p = setupGeneric(mAdult);
p = parametersChemostat(p);
p.L = 100;
p.tEnd=1000; %days

%make sure there is no default quadratic or declining HTL
setHTL(0,0,false,false);
%unicellulars already receive HTL mortality from multicellulars
phyto=zeros(1,10);
%multicellular mort definition
if options==1
    zoomort=repmat(mortHTL,length(p.ixStart(2):p.ixEnd(end)),1)';
    %merge unicellulars and multicellulars
    HTLmort=[phyto,zoomort];
elseif options==2
    mass=p.m(p.ixStart(2):end);
    pHTL=(1./(1+(mass./mortm).^(-2)));
    zoomort=mortHTL.*pHTL;
    HTLmort=[phyto,zoomort];
end

%call custom library
calllib(loadNUMmodelLibrary(), 'f_setmorthtl', ...
    double(HTLmort) );

sim=simulateChemostat(p, p.L); 

end