%
% Returns the trophic level of each cell groups in a vector lambda
%
% In:
%  p - a structure with parameters from a setup
%  rates - a structure with jDOC, jLreal and jF rates (from getRates)
%  lat,lon - coordinates (only for global simulation)
%
% Out:
%  lambda - the trophic level vector of all size groups that are not
%          nutrients.
%

function [lambda]=calcTrophicLevel(sim,rates,lat,lon)

arguments
    sim struct
    rates struct
    lat double = []; %only for global simulation
    lon double = []; %only for global simulation
end

%Get the size preference matrix of all size groups that are not nutrients.
theta = getTheta(sim.p);


%cell rates
jPP=rates.jDOC+rates.jLreal; %Uptake of Primary Production 
jF=rates.jFreal; %Uptake of Food

%Biomass
switch sim.p.nameModel
    case 'chemostat'
        B=sim.B;
    case 'watercolumn'
        B=squeeze(sum(sim.B(:,:,:),2));
    case 'global'
        idx = calcGlobalWatercolumn(lat,lon,sim);
        B=squeeze(sum(sim.B(:, idx.x, idx.y,sim.bathy(idx.x,idx.y,:)==1,:),4));
end

%
%Calculate the trophic level
%

%Sort the masses in ascending order
[~,mSorted]=sort(sim.p.m(sim.p.idxB:end));

lambda=B(1,:)*0;
lambda(mSorted(1))=1;
for i=2:length(lambda)
    idx=mSorted(i);
    if sum(theta(idx,:))==0
        lambda(idx)=1; %for cells not eating
    else
        lambda(idx)=(jPP(idx)+jF(idx)*(1+sum(theta(idx,:).*lambda.*B(end,:))/sum(theta(idx,:).*B(end,:))))/(jPP(idx)+jF(idx)); %trophic level
    end
end 



