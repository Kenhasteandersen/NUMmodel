%
% Simulate the POM in a Watercolumn   ---   Author: Cécile Decker
% In:
%  sim - The result of a simulation
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
    options.lon = -10;
    options.Z = 10^3*[0.0050 0.0150 0.0275 0.0450 0.0650 0.0875 0.1175 0.1600 0.2225 0.3100 0.4350 0.6100 0.8475 1.1600 1.5425];
end

rho_DOC = 1/3;
rho_DIC = 1-rho_DOC;
Gamma = 0.07; % remin rate (1/day) (Serra-Pompei 2022)
rho_CN = 5.68;
depthProductiveLayer = 50;

%Loading the 3 vectors (distrib, sinking velocity and radius, from Anton's
%studies (2D - density and size dependant POM) --> possible to change these ones easily
L = load('POM\data.mat');
POM.w = sum(L.M .* L.w)./sum(L.M);
%POM.w = sum(L.M .* L.w)/sum(sum(L.M)); %strange sinking velocity distribution ...
POM.r = L.r;
POM.distrib = sum(L.M .* L.w)/sum(sum(L.M .* L.w));
POM.m = sum(L.M);%/sum(sum(L.M));


%Calulation of export flux of POM from the productive layer

switch sim.p.nameModel
     
    case 'chemostat'

        %i = sim.p.ixPOM;
        % Part for generalists POM
        iStart = sim.p.ixStart(1) - sim.p.idxB +1;
        iEnd = sim.p.ixEnd(1) - sim.p.idxB +1;
        %B_0 = sum(sim.B(:,sim.p.ixStart(i):sim.p.ixEnd(i)-sim.p.idxB+1),2);
        % need to sum the biomass of POM of NUMmodel of each size group and multiply by the distrib export flux ...
        %B_i0 = (B_0*POM.distrib).';
        B_i0 = (sim.B(:,iStart:iEnd)*sim.rates.jPOM(iStart:iEnd)*POM.distrib)./POM.w;
        
        Z = options.Z(find(options.Z>depthProductiveLayer));
        POM.B_i = zeros(length(sim.t),length(Z),length(POM.distrib));
        POM.JDOC = zeros(length(sim.t),length(Z));
        POM.JDIC = zeros(length(sim.t),length(Z));
        POM.JN = zeros(length(sim.t),length(Z));
        %calculation for each depth
        for idx_z=1:length(Z)
            z = Z(idx_z);
            A = ones(length(POM.w),1)./POM.w.';
            POM.B_i(:,idx_z,:) = exp(-(Gamma*z*A)).'.*B_i0;
        
            %  Calculate the flux of Carbon and Nitrogen given from POM through depth
        
            POM.JDOC(:,idx_z) = rho_DOC*Gamma*squeeze(sum(POM.B_i(:,idx_z,:),3));
            POM.JDIC(:,idx_z) = rho_DIC*Gamma*squeeze(sum(POM.B_i(:,idx_z,:),3));
            POM.JN(:,idx_z) = (1/rho_CN)*Gamma*squeeze(sum(POM.B_i(:,idx_z,:),3));
        
        end
        POM.Z = Z;

    case 'watercolumn'

        %layer = 1; %for the moment we basicaly take the top layer
        %time = 150; %for the moment we simply take the 150th day of the simulation
        %i = sim.p.ixPOM;
        % Part for generalists POM
        iStart = sim.p.ixStart(1) - sim.p.idxB +1;
        iEnd = sim.p.ixEnd(1) - sim.p.idxB +1;
        iZ = (find(sim.z<=depthProductiveLayer));
        %B_0 = squeeze(sum(sim.B(layer,sim.p.ixStart(i):sim.p.ixEnd(i)-sim.p.idxB+1,:)));
        % need to sum the biomass of POM of NUMmodel of each size group and multiply by the distrib export flux ...
        %B_i0 = zeros(length(POM.distrib), length(B_0));
        %B_i0 = (B_0*POM.distrib).';
        B_i0 = squeeze(sum(sum((sim.B(iZ,iStart:iEnd,:).*sim.rates.jPOM(iZ,iStart:iEnd,:)))))*POM.distrib./POM.w;

        Z = sim.z(find(options.Z>depthProductiveLayer));
        POM.B_i = zeros(length(sim.t),length(Z),length(POM.distrib));
        POM.JDOC = zeros(length(sim.t),length(Z));
        POM.JDIC = zeros(length(sim.t),length(Z));
        POM.JN = zeros(length(sim.t),length(Z));
        %calculation for each depth
        for idx_z=1:length(Z)
            z = sim.z(idx_z);
            A = ones(length(POM.w),1)./POM.w.';
            POM.B_i(:,idx_z,:) = exp(-(Gamma*z*A)).'.*B_i0;
        
            %  Calculate the flux of Carbon and Nitrogen given from POM
            %  through depth

            POM.JDOC(:,idx_z) = rho_DOC*Gamma*squeeze(sum(POM.B_i(:,idx_z,:),3));
            POM.JDIC(:,idx_z) = rho_DIC*Gamma*squeeze(sum(POM.B_i(:,idx_z,:),3));
            POM.JN(:,idx_z) = (1/rho_CN)*Gamma*squeeze(sum(POM.B_i(:,idx_z,:),3));
        
        end
        POM.Z = Z;

    case 'global'

        %layer = 1; %for the moment we basicaly take the top layer
        %time = 150; %for the moment we simply take the 150th day of the simulation
        %i = sim.p.ixPOM;
        % Part for generalists POM
        iStart = sim.p.ixStart(1) - sim.p.idxB +1;
        iEnd = sim.p.ixEnd(1) - sim.p.idxB +1;
        iZ = (find(sim.z<=depthProductiveLayer));
        %epsilon = 2;
        lat = find(sim.x>options.lat,1);
        lon = find(sim.y>options.lon,1);
        %B_0 = squeeze(sum(sim.B(lat,lon,layer,sim.p.ixStart(i):sim.p.ixEnd(i)-sim.p.idxB+1,:)));
        % need to sum the biomass of POM of NUMmodel of each size group and multiply by the distrib export flux ...
        %B_i0 = zeros(length(POM.distrib), length(B_0));
        %B_i0 = (B_0*POM.distrib).';
        %B_i0 = zeros(length(sim.t), length(sim.x), length(sim.y), length(POM.distrib));
        Z = sim.z(find(options.Z>depthProductiveLayer));
        POM.B_i = zeros(length(sim.t), length(sim.x), length(sim.y), length(Z),length(POM.distrib));
        POM.JDOC = zeros(length(sim.t), length(sim.x), length(sim.y), length(Z));
        POM.JDIC = zeros(length(sim.t), length(sim.x), length(sim.y), length(Z));
        POM.JN = zeros(length(sim.t), length(sim.x), length(sim.y), length(Z));

        % Loop through lat and lon
        for ix=1:length(sim.x)
            for iy=1:length(sim.y)
                B_i0 = squeeze(sum(sum((sim.B(ix,iy,iZ,iStart:iEnd,:).*sim.rates.jPOM(ix,iy,iZ,iStart:iEnd,:)),3),4))*POM.distrib./POM.w;
        
                %calculation for each depth
                for idx_z=1:length(Z)
                    z = sim.z(idx_z);
                    A = ones(length(POM.w),1)./POM.w.';
                    POM.B_i(:,ix,iy,idx_z,:) = exp(-(Gamma*z*A)).'.*B_i0;
                
                    %  Calculate the flux of Carbon and Nitrogen given from POM
                    %  through depth
        
                    POM.JDOC(:,ix,iy,idx_z) = rho_DOC*Gamma*squeeze(sum(POM.B_i(:,idx_z,:),3));
                    POM.JDIC(:,ix,iy,idx_z) = rho_DIC*Gamma*squeeze(sum(POM.B_i(:,idx_z,:),3));
                    POM.JN(:,ix,iy,idx_z) = (1/rho_CN)*Gamma*squeeze(sum(POM.B_i(:,idx_z,:),3));
                
                end
            end
        end

        POM.Z = Z;


    otherwise
        error('Simulation type %s not supported.', sim.p.nameModel);

end

%% The END

% Maybe adding some basic plots
sim.POM = POM;
