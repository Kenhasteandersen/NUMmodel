%
% Simulate the POM in a Watercolumn   ---   Author: Cécile Decker
% In:
%  sim - The result of a Watercolumn simulation
%
% Options:
%  options - for the global simulation, we need to choose a certain lat and
%            lon where we want to see the POM flux. It can be changed after
%            to be seen worldwide, bot one has to be careful on the tricky
%            dimensions.
%
% Out:
%  sim - simulation object with POM added but also matching sinking
%        velocity and mass
%
function sim = calcPOM(sim,options)

arguments
    sim = simulateWatercolumn;
    options.lat = 60;
    options.lon = -20;
end

rho_DOC = 1/3;
rho_DIC = 1-rho_DOC;
Gamma = 0.07; % remin rate (1/day) (Serra-Pompei 2022)
rho_CN = 5.68;

%Loading the 3 vectors (distrib, sinking velocity and radius, from Anton's
%studies (2D - density and size dependant POM) --> possible to change these ones easily
L = load('POM\data.mat');
POM.w = sum(L.M .* L.w)/sum(sum(L.M));
POM.r = L.r;
POM.distrib = sum(L.M .* L.w)/sum(sum(L.M .* L.w));



%Calulation of total mass of POM at the surface layer

switch sim.p.nameModel
    
    case 'chemostat'

        i = sim.p.ixPOM;
        B_0 = sum(sim.B(:,sim.p.ixStart(i):sim.p.ixEnd(i)-sim.p.idxB+1),2);
        % need to sum the biomass of POM of NUMmodel of each size group and multiply by the distrib export flux ...
        B_i0 = (B_0*POM.distrib).';
        
        Z = load('POM\Z.mat').Z;
        POM.B_i = zeros(length(Z),length(POM.distrib),length(sim.t));
        POM.JDOC = zeros(length(Z),length(sim.t));
        POM.JDIC = zeros(length(Z),length(sim.t));
        POM.JN = zeros(length(Z),length(sim.t));
        %calculation for each depth
        for idx_z=1:length(Z)
            z = Z(idx_z);
            POM.B_i(idx_z,:,:) = exp(-(POM.w/Gamma)*z).'.*B_i0;
        
            %  Calculate the flux of Carbon and Nitrogen given from POM through depth
        
            POM.JDOC(idx_z,:) = rho_DOC*Gamma*squeeze(sum(POM.B_i(idx_z,:,:)));
            POM.JDIC(idx_z,:) = rho_DIC*Gamma*squeeze(sum(POM.B_i(idx_z,:,:)));
            POM.JN(idx_z,:) = (1/rho_CN)*Gamma*squeeze(sum(POM.B_i(idx_z,:,:)));
        
        end

    case 'watercolumn'

        layer = 1; %for the moment we basicaly take the top layer
        %time = 150; %for the moment we simply take the 150th day of the simulation
        i = sim.p.ixPOM;
        B_0 = squeeze(sum(sim.B(layer,sim.p.ixStart(i):sim.p.ixEnd(i)-sim.p.idxB+1,:)));
        % need to sum the biomass of POM of NUMmodel of each size group and multiply by the distrib export flux ...
        B_i0 = zeros(length(POM.distrib), length(B_0));
        B_i0 = (B_0*POM.distrib).';
        
        POM.B_i = zeros(length(sim.z),length(POM.distrib),length(sim.t));
        POM.JDOC = zeros(length(sim.z),length(sim.t));
        POM.JDIC = zeros(length(sim.z),length(sim.t));
        POM.JN = zeros(length(sim.z),length(sim.t));
        %calculation for each depth
        for idx_z=1:length(sim.z)
            z = sim.z(idx_z);
            POM.B_i(idx_z,:,:) = (exp(-(POM.w/Gamma)*z).'.*B_i0);
        
            %  Calculate the flux of Carbon and Nitrogen given from POM through depth
        
            POM.JDOC(idx_z,:) = rho_DOC*Gamma*squeeze(sum(POM.B_i(idx_z,:,:)));
            POM.JDIC(idx_z,:) = rho_DIC*Gamma*squeeze(sum(POM.B_i(idx_z,:,:)));
            POM.JN(idx_z,:) = (1/rho_CN)*Gamma*squeeze(sum(POM.B_i(idx_z,:,:)));
        end

    case 'global'

        layer = 1; %for the moment we basicaly take the top layer
        %time = 150; %for the moment we simply take the 150th day of the simulation
        i = sim.p.ixPOM;
        epsilon = 2;
        lat = find(abs(sim.x-options.lat)<epsilon,1);
        lon = find(abs(sim.y-options.lon)<epsilon,1);
        B_0 = squeeze(sum(sim.B(lat,lon,layer,sim.p.ixStart(i):sim.p.ixEnd(i)-sim.p.idxB+1,:)));
        % need to sum the biomass of POM of NUMmodel of each size group and multiply by the distrib export flux ...
        B_i0 = zeros(length(POM.distrib), length(B_0));
        B_i0 = (B_0*POM.distrib).';
        
        POM.B_i = zeros(length(sim.z),length(POM.distrib),length(sim.t));
        POM.JDOC = zeros(length(sim.z),length(sim.t));
        POM.JDIC = zeros(length(sim.z),length(sim.t));
        POM.JN = zeros(length(sim.z),length(sim.t));
        %calculation for each depth
        for idx_z=1:length(sim.z)
            z = sim.z(idx_z);
            POM.B_i(idx_z,:,:) = (exp(-(POM.w/Gamma)*z).'.*B_i0);
        
            %  Calculate the flux of Carbon and Nitrogen given from POM through depth
        
            POM.JDOC(idx_z,:) = rho_DOC*Gamma*squeeze(sum(POM.B_i(idx_z,:,:)));
            POM.JDIC(idx_z,:) = rho_DIC*Gamma*squeeze(sum(POM.B_i(idx_z,:,:)));
            POM.JN(idx_z,:) = (1/rho_CN)*Gamma*squeeze(sum(POM.B_i(idx_z,:,:)));
        end


    otherwise
        error('Simulation type %s not supported.', sim.p.nameModel);

end

%% The END

% Maybe adding some basic plots
sim.POM = POM;
