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

function [lambda]=calcTrophicLevel(p,rates)

arguments
    p struct
    rates struct
end

%Get the size preference matrix of all size groups that are not nutrients.
theta = getTheta(p);

%cell rates
jPP=rates.jDOC+rates.jLreal; %Uptake of Primary Production 
jF=rates.jF; %Uptake of Food

%
%Calculate the trophic level
%
lambda=p.u0(p.idxB:end)*0;
lambda(1)=1;
for i=2:length(lambda)
    if sum(theta(i,1:i-1))==0
        lambda(i)=1; %for cells not eating
    else
        lambda(i)=(jPP(i)+jF(i)*(1+sum(theta(i,1:i-1).*lambda(1:i-1)))/sum(theta(i,1:i-1)))/(jPP(i)+jF(i)); %trophic level
    end
end 



