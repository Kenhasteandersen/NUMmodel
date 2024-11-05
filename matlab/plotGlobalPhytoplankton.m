%
% Make a plot of bacteria, phytoplankton, zooplankton, and the diatom ratio.
%
% Requires that the library is loaded with a setup which matches sim and
% that it is setup with bParallel=true.
% 
function sim = plotGlobalPhytoplankton(sim, options)
arguments
    sim struct;
    options.sProjection = 'fast'; %projection to use. Defaults to 'fast'. Other projections
%               requires that the mapping toolbox is installed. 
%               Good projection is 'eckert4'.
end
sLibName = loadNUMmodelLibrary();
ixTime = find(sim.t>(max(sim.t)-365)); %nTime = length(sim.t(sim.t >= max(sim.t-365))); % Just do the last year
% Get grid volumes:
%load(sim.p.pathGrid,'dv','dz','dx','dy');
dz = sim.dznom;
ix = ~isnan(sim.N(1,:,:,1)); % Find all relevant grid cells


Bphyto = zeros(length(sim.t), length(sim.x), length(sim.y));
Bbacteria = Bphyto;
Bzoo = Bphyto;
BphytoDiatoms = Bphyto;
BphytoOthers = Bphyto;

nX = length(sim.x);
nY = length(sim.y);
nZ = length(sim.z);
%
% Extract fields from sim:
%
N = sim.N;
DOC = sim.DOC;
if isfield(sim.p,'idxSi')
    Si = sim.Si;
else
    Si = 0;
end
B = sim.B;
L = sim.L;
T = sim.T;
p = sim.p;
%
% find indices of diatoms and generalists:
%
ixGeneralists = p.typeGroups==1 | p.typeGroups==5;
ixDiatoms = p.typeGroups==3 | p.typeGroups==4;

parfor iTime = ixTime
    for i = 1:nX
        for j = 1:nY
            for k = 1:nZ
                if ~isnan(N(iTime,i,j,k))
                    if isfield(p,'idxSi')
                        u = [squeeze(N(iTime,i,j,k)), ...
                            squeeze(DOC(iTime,i,j,k)), ...
                            squeeze(Si(iTime,i,j,k)), ...
                            squeeze(B(iTime,i,j,k,:))'];
                    else
                        u = [squeeze(N(iTime,i,j,k)), ...
                            squeeze(DOC(iTime,i,j,k)), ...
                            squeeze(B(iTime,i,j,k,:))'];
                    end
                    [Bphytotmp, Bzootmp, Bbacteriatmp] = ...
                        calcPhytoZoo(p, u, L(iTime,i,j,k), T(iTime,i,j,k), sLibName);
                    
                    conv = dz(k);%squeeze(dz(i,j,k));
                    Bphyto(iTime,i,j) = Bphyto(iTime,i,j) + sum(Bphytotmp)*conv; % mgC/m2
                    Bbacteria(iTime,i,j) = Bbacteria(iTime,i,j) + sum(Bbacteriatmp)*conv; % mgC/m2
                    Bzoo(iTime,i,j) = Bzoo(iTime,i,j) + sum(Bzootmp)*conv; % mgC/m2

                    BphytoDiatoms(iTime,i,j) = BphytoDiatoms(iTime,i,j)+ sum(Bphytotmp(ixDiatoms))*conv;
                    BphytoOthers(iTime,i,j) = BphytoOthers(iTime,i,j) + sum(Bphytotmp(ixGeneralists))*conv;
                end
            end
        end
    end
end
sim.Bphyto = Bphyto;
sim.Bbacteria = Bbacteria;
sim.Bzoo = Bzoo;
sim.Diatomratio = BphytoDiatoms ./ (BphytoDiatoms+BphytoOthers);
sim.Bdiatoms = BphytoDiatoms;
%%
% Make plot:
%
clf
tiledlayout(4,1)

nexttile
panelGlobal(sim.x, sim.y, mean(sim.Bbacteria(ixTime,:,:),1), sTitle="Bacteria", sUnits="mg_C/m^2", sProjection=options.sProjection)


nexttile
panelGlobal(sim.x, sim.y, mean(sim.Bphyto(ixTime,:,:),1), sTitle="Phytoplankton", sUnits="mg_C/m^2", sProjection=options.sProjection)

nexttile
panelGlobal(sim.x, sim.y, mean(sim.Bzoo(ixTime,:,:),1), sTitle="Zooplankton", sUnits="mg_C/m^2", sProjection=options.sProjection)

nexttile
panelGlobal(sim.x, sim.y, mean(sim.Diatomratio(ixTime,:,:),1), sTitle="Diatom ratio", sUnits="", sProjection=options.sProjection)
clim([0 1])
