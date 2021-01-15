function sim = calcGlobalFunction(sim)
% Get grid volumes:
load(sim.p.pathGrid,'dv','dx','dy','dz');
ix = ~isnan(sim.N(:,:,:,1)); % Find all relevant grid cells

%
% Primary production:
%
if ~isfield(sim,'Cnet')
    for i = 1:length(sim.t)
        sim.Cnet(:,:,:,i) = calcGlobalCnet(sim,i); % mugC/l/day
    end
end
%
% Global totals
%
    function tot = calcTotal(u)
        tot = sum(u(ix(:)).*dv(ix(:)))*10*10; % mug/day
    end

for i = 1:length(sim.t)
    sim.Ntotal(i) = calcTotal(sim.N(:,:,:,i));
    sim.DOCtotal(i) = calcTotal(sim.DOC(:,:,:,i)); % mugC
    sim.Btotal(i) = 0;
    for j = 1:sim.p.n
        sim.Btotal(i) = sim.Btotal(i) + calcTotal(sim.B(:,:,:,j,i)); % mugC
    end
    sim.CnetTotal(i) = calcTotal(sim.Cnet(:,:,:,i)); % mugC/day
end
%
% Watercolumn totals:
%
for i = 1:length(sim.x)
    for j = 1:length(sim.y)
        for k = 1:length(sim.t)       
            % gC per m2/year:
            sim.CnetPerArea(i,j,k) = sum(squeeze(sim.Cnet(i,j,:,k)).*squeeze(dz(i,j,:)))*365*10*10*1000; 
        end
    end
end
%
% Annual global totals:
%
sim.CnetPerAreaAnnual = zeros(length(sim.x), length(sim.y), floor(sim.t(end)/365));
for i = 1:sim.t(end)/365
    ixTime = sim.t>365*(i-1) & sim.t<=365*i;
    sim.CnetTotalAnnual(i) = sum(sim.CnetTotal(ixTime))/length(ixTime); % mugC/l/day
    sim.CnetPerAreaAnnual(:,:,i) = sum(sim.CnetPerArea(:,:,ixTime),3)/length(ixTime); % gC/m2/yr
end
sim.CnetTotalAnnual = sim.CnetTotalAnnual*365/1000/1000/1000/1000; % GTon carbon/year

end