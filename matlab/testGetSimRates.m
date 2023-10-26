lat_to_find=55; lon_to_find=-40;
noYears=20;
final_years=3;
mHTL = 1;
mortHTL = .15;
bHTLdecline = true;
bHTLquadratic = true;   
runningTime = noYears*365;
newSinkingPOM = 0.78622;
newRemin2= logspace(log10(0.01), log10(0.5), 10);
newParameter=newRemin2;
load(['sim_remin2_',num2str(newRemin2(1)),'lat',num2str(lat_to_find),'.mat'])
sim.p = parametersWatercolumn( setupNUMmodel );
        setHTL(mortHTL, mHTL, bHTLquadratic, bHTLdecline)
        setSinkingPOM(sim.p, newSinkingPOM)
        sim.p.tEnd = runningTime;

% simR=getSimRates(sim);
timePeriod=(noYears-final_years)*365+1:runningTime;
% Bphyto=zeros(length((noYears-final_years)*365+1:runningTime),20,length(sim.z));
% 
% Bph_orig=calcPhyto(sim,simR);
% Bph=Bph_orig;
z = sim.z + 0.5*sim.dznom;
% Make a layer at z = 0 with the same value as in the first grid point:
%
t = sim.t;
z = [0; sim.z];
%%
Bphyto_avg=zeros(noYears,12);
Bph_param=zeros(length(newParameter),12);
for iparam=1:length(newParameter)
    mat_to_load=['sim_remin2_',num2str(newRemin2(iparam)),'lat',num2str(lat_to_find),'.mat'];
    load(mat_to_load)
    simR=getSimRates(sim);
    Bph_orig=calcPhyto(sim,simR);
    Bph=Bph_orig;

    Bph(:,2:length(z),:) = Bph;
    Bph(:,1,:) = Bph(:,2,:);
    Bph_allSizes=squeeze(sum(Bph(timePeriod,:,:),3));
    for i= (noYears-final_years)+1:noYears
        Bphyto_avg(i,:)=reshapeCellToArrayAvg(Bph_allSizes,i);
    end

    Bph_param(iparam,:)=squeeze(mean(Bphyto_avg(noYears-3+1:noYears,:),1));
end
save(['Bph_param_remin2_lat',num2str(lat_to_find),'.mat'],'Bph_param');    
