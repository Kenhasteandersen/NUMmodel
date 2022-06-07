%
% Simulate the POM in a Watercolumn   ---   Author: Cécile Decker
% In:
%  sim - The result of a Watercolumn simulation
%
% Options:
%
% Out:
%  sim - simulation object with POM added but also matching sinking
%        velocity and mass
%
function sim = calcPOM(sim)

arguments
    sim = simulateWatercolumn;
end

rho_DOC = 1/3;
rho_DIC = 1-rho_DOC;
Gamma = 0.07; % remin rate (1/day) (Serra-Pompei 2022)
rho_CN = 5.68;

%% everything concerning the POM spectrum in the Watercolumn

%Loading the 3 vectors (distrib, sinking velocity and radius, from Anton's
%studies (2D - density and size dependant POM) --> possible to change these ones easily
L = load('POM_data\export_cecile.mat');
POM.w = L.w;
POM.r = L.r;
POM.distrib = L.export_frac;

%plot(log(POM.r), POM.distrib) % plot(POM.r, POM.distrib)
%title("Distribution of POM in the surface layer, calculated with Anton's 2D representation")
%ylabel('distribution')
%xlabel('log10 of sizes (\mu m)') %xlabel('sizes (\mu m)')

%plot(POM.r, POM.w)
%xlabel('sizes (\mu m)')
%title("Sinking velocity of POM calculated with Anton's 2D representation (mass weighted)")
%ylabel('sinking velocity (m/d)')

%Calulation of total mass of POM at the surface layer
% Need to be rebuilt to match other types of simulation (Chemostat/General)
layer = 1; %for the moment we basicaly take the top layer
time = 150; %for the moment we simply take the 150th day of the simulation
%we take the index correspondant to the POM : sim.p.ixStart(sim.p.ixPOM):sim.p.ixEnd(sim.p.ixPOM)
%B_0 = sum(sim.B(layer,61:end,time));
i = sim.p.ixPOM;
B_0 = sum(sim.B(layer,sim.p.ixStart(i):sim.p.ixEnd(i)-sim.p.idxB+1,time));
% need to sum the biomass of POM of NUMmodel of each size group and multiply by the distrib export flux ...
B_i0 = B_0*POM.distrib;

%calculation for each depth
for idx_z=1:length(sim.z)
    z = sim.z(idx_z);
    %B_i0
    %POM.w
    B_i(idx_z,:) = B_i0.*exp(-(POM.w/Gamma)*z);

    %%  Calculate the flux of Carbon and Nitrogen given from POM through depth

    POM.JDOC(idx_z) = rho_DOC*Gamma*sum(B_i(idx_z,:));
    POM.JDIC(idx_z) = rho_DIC*Gamma*sum(B_i(idx_z,:));
    POM.JN(idx_z) = (1/rho_CN)*Gamma*sum(B_i(idx_z,:));

end

%% The END

% Maybe adding some basic plots
sim.POM = POM;

% Plots
%figure(1)
%clf
%nexttile
%plot(POM.JDOC', -sim.z,'linewidth', 2)%, '-')
%title('flux of DOC, from the POM')
%xlabel('DOC flux')
%ylabel('depth')
%nexttile
%plot(POM.JDIC, -sim.z)
%title('flux of DIC, from the POM')
%xlabel('DIC flux')
%ylabel('depth')
%nexttile
%plot(POM.JN, sim.z)
%title('flux of Nitrogen, from the POM')
%xlabel('Nitrogen flux')
%ylabel('depth')

%figure(2)
%clf
%nexttile
