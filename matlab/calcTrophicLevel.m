%
% Returns the trophic level of each cell groups in a vector lambda
%
% In:
%  p - a structure with parameters from a setup
%  rates - a structure with jDOC, jLreal and jF rates (from getRates)
%
% Out:
%  lambda - the trophic level vector of all size groups that are not
%          nutrients.
%

function [lambda]=calcTrophicLevel(sim,rates)

arguments
    sim struct
    rates struct
end

%Get the size preference matrix of all size groups that are not nutrients.
theta = getTheta(sim.p);


%cell rates
jPP=rates.jDOC+rates.jLreal; %Uptake of Primary Production 
jF=rates.jFreal; %Uptake of Food

%
%Calculate the trophic level
%

%Sort the masses in ascending order
[~,mSorted]=sort(sim.p.m(sim.p.idxB:end));

lambda=sim.B(1,:)*0;
lambda(mSorted(1))=1;
for i=2:length(lambda)
    idx=mSorted(i);
    if sum(theta(idx,:))==0
        lambda(idx)=1; %for cells not eating
    else
        lambda(idx)=(jPP(idx)+jF(idx)*(1+sum(theta(idx,:).*lambda.*sim.B(end,:))/sum(theta(idx,:).*sim.B(end,:))))/(jPP(idx)+jF(idx)); %trophic level
    end
end 



