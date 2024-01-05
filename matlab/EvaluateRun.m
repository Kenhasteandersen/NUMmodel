function status = EvaluateRun(sim)

sim.Bpico = sim.BpicoAnnualMean;
r = calcRadiusGroups(sim.p);

figure(1)
clf
tiles = tiledlayout(2,4);
tiles.TileIndexing = 'columnmajor';
tiles.TileSpacing = 'tight';
ylimLatitude = [-65 55];

% Pico data:
fitPico = PicophytoFiguresAtlantic(sim,r);
fprintf('Pico fit:    %f\n', fitPico);
% POPcorn data:
fitPopcorn = GOPOPCORNfigures(sim,r);
fprintf('Popcorn fit: %f\n', fitPopcorn);
% Copepod data:
fitCopepods = AMTcopepods(sim);
fprintf("Copepods fit: %f\n", fitCopepods);
%NPP:
%figure(2)
NPPcafeVsNUMfinal(sim)
NPPEutro = printStats(0,-119,sim);
NPPOligo = printStats(22,158,sim);
NPPSeasonal = printStats(60,-15,sim);

status = [fitPopcorn, fitPico, fitCopepods, NPPEutro, NPPOligo, NPPSeasonal];

    function NPP = printStats(lat,lon,sim)
        ix = calcGlobalWatercolumn(lat,lon,sim);
        ixTime = sim.t>max(sim.t)-365;
        %fprintf('Lat %i, lon %i:\n',lat,lon);
        NPP = mean(spatialaverage(ix,sim.ProdNetAnnual(1,:,:)));
        fprintf('Prod net %5.0f mgC/m2/day\n', NPP);
        %fprintf('Pico      %2.2f mgC/m2\n', mean(spatialaverage(ix,sim.BpicoAnnualMean(1, :,:))));
        r = calcRadius(sim.p.m(sim.p.idxB:end));
        %POC = mean(sum(sim.B(sim.t>sim.t(end)-365,ix.x,ix.y,1,r>0.7 & r<30),5),1);
        %fprintf('POC        %2.2f mgC/m3\n', POC);
        %fprintf('----\n');
    end

    function phi = spatialaverage(ix,field)
        phi = (...
            field(:,ix.x, ix.y) + ...
            field(:,ix.x-1, ix.y) + ...
            field(:,ix.x+1, ix.y) + ...
            field(:,ix.x, ix.y+1) + ...
            field(:,ix.x-1, ix.y+1) + ...
            field(:,ix.x+1, ix.y+1) + ...
            field(:,ix.x, ix.y-1) + ...
            field(:,ix.x-1, ix.y-1) + ...
            field(:,ix.x+1, ix.y-1) )/9;
    end

% Converts generalists and diatoms mass to radius in mum
% Converts copepod body-mass to prosome length/2 in mum, based on
% prosome length(um) to body-mass relationship for copepods
% from Chisholm and Roff (1990)
% -------treat POM as generalists------- MUST BE ADDRESSED
    function r = calcRadiusGroups(p)

        rho = 0.4*1e6*1e-12; % mug/cm3 (Andersen et al 2016; rho = m/V = 0.3*d^3/(4/3*pi*(d/2)^3) )

        for iGroup = 1:p.nGroups
            ix = (p.ixStart(iGroup):p.ixEnd(iGroup))-p.idxB+1;
            m(ix) = p.m(ix+p.idxB-1);
            % ix = (p.ixStart(iGroup):p.ixEnd(iGroup));
            % m = p.m(ix);

            %
            % Generalists:
            %
            if ((p.typeGroups(iGroup)==1) || (p.typeGroups(iGroup)==5)|| (p.typeGroups(iGroup)==100))
                r(ix) = (3/(4*pi)*m(ix)/rho).^(1/3);
            end
            %
            % Diatoms:
            %
            if ((p.typeGroups(iGroup)==3) || (p.typeGroups(iGroup)==4))
                v = 0.6;
                r(ix) = (3/(4*pi)*m(ix)/rho/(1-v)).^(1/3);
            end
            if ((p.typeGroups(iGroup)==10) || (p.typeGroups(iGroup)==11))
                % prosome length divided by 2 to be equivalent to radius
                r(ix)= ( m(ix)/(0.73*0.48*exp(-16.41))^(1/2.74) )/2;
            end

        end
    end

% interpolateToNUM and recenter
    function npp_interp=interpolateToNUM(sim,npp_annual_avg)
        lat=zeros(1080,1);
        lat(1)=1/6;
        for i=2:1080
            lat(i)=lat(i-1)+1/6;
        end
        lat=lat-90;
        %
        lon=zeros(2160,1);
        lon(1)=1/6;
        for i=2:2160
            lon(i)=lon(i-1)+1/6;
        end

        % create meshgrid with CAFE's coordinates
        [lon1,lat1]=meshgrid(lon,lat);

        %create meshgrid with NUM's coordinates
        [lon2,lat2]=meshgrid(sim.x,sim.y);

        % interpolate NPP CAFE dimensions to NUM's
        nppCAFEinterp=interp2(lon1,lat1,squeeze(npp_annual_avg),lon2,lat2);

        orig0=(sim.ProdNetAnnual(end,:,:)); %NUM NPP
        orig=squeeze(orig0)'; % flip NUM npp data

        nppCAFEinterpRecenter= [nppCAFEinterp(:,65:end),nppCAFEinterp(:,1:64)];
        % make land white
        nppCAFEinterpRecenter(nppCAFEinterpRecenter==-9999)=NaN;
        orig(orig==0)=NaN;


        npp_interp=nppCAFEinterpRecenter;
    end



%--------------------------------------------------------------------------
% Plots of GO-POPCORN data and comparison with model output (sim.Bmicro)
% data: https://doi.org/10.1038/s41597-022-01809-1
% this is only the code, the model results here don't fit
% (the biomass is too low, but it is probably better in the newest version)
% copepods and detritus should also be include in Bmicro-POC_avg_uM (?)
%
%--------------------------------------------------------------------------
%                          DATA PROCESSING
%-   -   - -   - -   - -   - -   - -   - -   - -   - - -   - - -   - - -
% *STEP 1*: Convert data table to array and convert units of POC to ugC/L
% *STEP 2*: Select size-classes of the model with r in [0.7/2, 30/2] um
% *STEP 3*: Remove the rows where POC is NaN
% *STEP 4*: Group data by cruise and subgroup by cruise_station
% *STEP 5*: Calculate the mean values at each station --> 'summary_Dep'
% *STEP 6*: Isolate only the last year of the model output and match the
%           model months to the observation months
% *STEP 7*: Find the model biomass at the same coordinates, depth and month
%           as the data, 'Pbiom_model'
% *STEP 8*: Calculate the mean value over depth at each coordinate,
%           'Pbiom_model_depth_avg'
% *STEP 9*: Apply the same grouping as the data (STEP 4) to the model
%           output from STEP 8 'summary_model'
% *STEP 10*: Plot transect
% *STEP 11*: For each cruise, plot observations and model output
% - -   - - -   - - -   - - -   - - -   - - -   - - -   - - -   - - -   - -
%--------------------------------------------------------------------------
% Arctic cruises are not included
    function fit = GOPOPCORNfigures(sim,r)


        load('dataGOPOPCORNver2.mat', 'dataGOPOPCORNver2')
        %%
        % sim=simFun1;

        % MODEL OUTPUT
        % load('simGlobalDevelopCompareFun.mat')
        % sim=simGlobalDevelopComapreFun;


        % convert data variables to vectors
        dataPOP=dataGOPOPCORNver2;
        cruise=table2array(dataPOP(:,"Cruise"));
        cruise_station=table2array(dataPOP(:,"Cruise_Station"));
        lat=table2array(dataPOP(:,"Latitude"));
        lat1=lat;
        lon=table2array(dataPOP(:,"Longitude"));
        lon1=lon;
        depth=table2array(dataPOP(:,"Depth"));
        month=table2array(dataPOP(:,"Month"));
        day=table2array(dataPOP(:,"Day"));
        year=table2array(dataPOP(:,"Year"));
        POCavg_uM=table2array(dataPOP(:,"POCavg_uM")).*12;  % converted to mugC/L

        %nutricline_1uM_Interp=table2array(dataPOP(:,18));

        lonWrapped = wrapTo180(sim.x); % convert longitude ordinates from [0,360[-->[-180,180]

        % The size range is 2.4ug-30um, so Biomass needs to be updated.
        % r = calcRadiusGroups(sim.p);
        % load('rGroups.mat')
        % r=rGroups;
        r_popcorn=find(r>=0.7/2 & r<=30/2);

        Bpopcorn=sim.B(:,:,:,:,r_popcorn); % mugC/L
        Bpopcorn_sum=squeeze(sum(Bpopcorn,5)); % sum all size-classes
        %%
        dataPOPc=[cruise cruise_station lat lon depth month POCavg_uM];
        nan_rows= any(isnan(dataPOPc(:,7)),2); %POCavg_uM=NaN
        % Remove rows where POC=NaN
        %...........................
        dataPOPclean = dataPOPc(~nan_rows, :);

        % if cruise_station==NaN then assign it to a new non-existing unique station
        for i=1: height(dataPOPclean)
            if isnan(dataPOPclean(i,2))
                dataPOPclean(i,2)=1000+i;
            end
        end

        % convert matrix to table
        dataPOP = array2table(dataPOPclean);
        % Create custom column names
        customColumnNames = {'Cruise', 'Cruise_Station', 'Latitude','Longitude','Depth','Month','POCavg_uM'};
        dataPOP.Properties.VariableNames = customColumnNames;

        %% -------------------------------------------------------------------
        % grouping based on cruise and subgrouping based on cruise station
        %-------------------------------------------------------------------

        grouping_cruise_cs=findgroups(dataPOP.Cruise, dataPOP.Cruise_Station);
        % create table same as dataPOP including the 'Grouping'
        temp_table=table(grouping_cruise_cs,dataPOP.Cruise,dataPOP.Latitude,dataPOP.Longitude,dataPOP.Depth,...
            dataPOP.Month,dataPOP.POCavg_uM,'VariableNames',{'Grouping','Cruise','Latitude','Longitude','Depth','Month','POCavg_uM'});
        %.................................................
        % average POC concentartion and depth at each subgroup.
        %.................................................
        summary_Dep=groupsummary(temp_table,'Grouping','mean'); % this gives us mean values at each station

        temp_array=table2array(temp_table);
        % remove arctic cruises: 1701,1709,1901,2018
        temp_array_noA=temp_array(1:2296,:); %2296
        myData=temp_array_noA;
        %%

        %---------------------------------------------------------------------
        % Find model Biomass for the same coordinates and depth with the data
        %---------------------------------------------------------------------
        lat=myData(:,3);
        lon=myData(:,4);
        depth=myData(:,5);
        month=myData(:,6);
        % Initialize Pbiom_model with zeros
        Pbiom_model = zeros(length(lat),1);
        Pbiom_model_depth_avg=Pbiom_model;

        for i=1:length(lat)

            [ ~, idx_lonG] = min( abs( lonWrapped-lon(i) ) );

            [ ~, idx_latG] = min( abs( sim.y-lat(i) ) );

            [ ~, idx_depth] = min( abs( sim.z-depth(i) ) );

            % Find the index of the corresponding month over the last year of the
            % simulation
            % [ ~, idx_mon] = min( abs( sim.t-(month(i)+144) ) );
            idx_mon = (length(sim.t)-12) + month(i);

            % Save the indices just in case
            ix_lon(i)=idx_lonG;
            ix_lat(i)=idx_latG;
            ix_dep(i)=idx_depth;
            ix_mon(i)=idx_mon;

            % Assign POC values for the indices above
            Pbiom_model(i)=Bpopcorn_sum(idx_mon(1),idx_lonG(1),idx_latG(1),idx_depth(1)); % in mugC/L (I think)
            Pbiom_model_depth_avg(i)= mean(Bpopcorn_sum(idx_mon(1),idx_lonG(1),idx_latG(1),idx_depth(1)),4)';
            % Pbiom_modelnew(i,:)=Bpopcorn_sum(:,idx_lonG(1),idx_latG(1),idx_depth); %(all month from the simulation)

        end

        depth_model=sim.z(ix_dep);
        lat_model=sim.y(ix_lat);
        lon_model=lonWrapped(ix_lon);
        month_model=mod(ix_mon,12)';
        for i=1:length(month_model)
            if month_model(i)==0
                month_model(i)=12;
            end
        end

        %%
        grouping_numbers=unique(myData(:,1));
        %(group_cruise_name(:,2)==cruise)
        sumPOC_group=zeros(1,length(grouping_numbers));
        sumPOC_model=zeros(1,length(grouping_numbers));

        count=zeros(1,length(grouping_numbers));
        count_model=zeros(1,length(grouping_numbers));

        depth_max=zeros(1,length(grouping_numbers));
        myData_wo_model=myData;
        myData_with_model=[myData lat_model lon_model depth_model Pbiom_model month_model];

        myData_sorted=sortrows(myData_with_model,1); % sorted data based on grouping
        myData_sorted_modelOnly=[myData_sorted(:,1) myData_sorted(:,8:11)];
        for ix_g=1:length(grouping_numbers)
            for j=1:length(myData_sorted)
                if myData_sorted(j,1)==ix_g
                    sumPOC_group(ix_g)=sumPOC_group(ix_g) +myData_sorted(j,7);
                    count(ix_g)=count(ix_g)+1;
                    lat_group(ix_g)=myData_sorted(j,3);
                    lon_group(ix_g)=myData_sorted(j,4);
                    cruise_name(ix_g)=myData_sorted(j,2);
                    month_group(ix_g)=myData_sorted(j,6);
                    % if depth_update~=myData_with_model(j,10)
                    %     sumPOC_model(ix_g)=sumPOC_model(ix_g) +myData_with_model(j,11);
                    %     count_model(ix_g)=count_model(ix_g)+1;
                    %     depth_update=myData_with_model(j,10);
                    % end
                    if myData(j,5)>depth_max(ix_g)
                        depth_max(ix_g)=myData(j,5); % max_depth at each grouping/station
                    end
                end
            end
        end


        myData_model_unique=unique(myData_sorted_modelOnly,"rows");
        myData_model_unique2=unique(myData_sorted_modelOnly(:,1:end-1),"rows");

        for ix_g=1:length(grouping_numbers)
            for j=1:length(myData_model_unique)
                if myData_model_unique(j,1)==ix_g
                    sumPOC_model(ix_g)=sumPOC_model(ix_g) +myData_model_unique(j,5);
                    count_model(ix_g)=count_model(ix_g)+1;
                    lat_model_new(ix_g)=myData_model_unique(j,2);
                    lon_model_new(ix_g)=myData_model_unique(j,3);
                    if myData_model_unique(j,4)>depth_max(ix_g)
                        depth_max_model(ix_g)=myData_model_unique(j,4); % max_depth at each grouping/station
                    end
                end
            end
        end



        for ix_g=1:length(grouping_numbers)
            meanPOC_group(ix_g)=sumPOC_group(ix_g)/count(ix_g);
            meanPOC_model(ix_g)=sumPOC_model(ix_g)/count_model(ix_g);
        end

        myData_final=[grouping_numbers cruise_name' lat_group' lon_group' month_group' meanPOC_group' lat_model_new' lon_model_new' meanPOC_model'];


        colNames = {'grouping','cruise name','lat','lon','month','mean POC','lat model','lon model','mean POC model'};
        finalDataTable = array2table(myData_final,'VariableNames',colNames);
        %%
        % cruise 46 has only Nan
        cruise_name=[7,9,13,18,28,1319,1418,1701,1709,1901,2018];

        % find indices of each cruise
        cruise1_ix=find(myData_final(:,2)==cruise_name(1));
        cruise2_ix=find(myData_final(:,2)==cruise_name(2));
        cruise3_ix=find(myData_final(:,2)==cruise_name(3));
        cruise4_ix=find(myData_final(:,2)==cruise_name(4));
        cruise5_ix=find(myData_final(:,2)==cruise_name(5));
        cruise6_ix=find(myData_final(:,2)==cruise_name(6));
        cruise7_ix=find(myData_final(:,2)==cruise_name(7));

        months=1:12;

        %%          TRANSECTS PLOT
        % figure
        % clf
        %    set(gcf,'color','w');
        %    surface(lon,lat,POCavg_uM,'EdgeColor','none')
        %    axis tight
        %    title('POC AVG uM')

        cmap=  [255, 0, 0; 255, 128, 0; 255, 255, 0;...
            128, 255, 0; 0, 255, 0; 0, 255, 128; 0, 255, 255;...
            0, 128, 255; 0, 0, 255; 128, 0, 255; 0, 210, 241;...
            255, 0, 255; 255, 0, 128]./255;

        latlim=[min(lat),max(lat)];
        lonlim=[min(lon),max(lon)];
        %axesm('miller','MapLatLimit',latlim,'MapLonLimit',lonlim,...
        %    'Frame','on','Grid','on','MeridianLabel','on','ParallelLabel','on','fontsize',8)
        %hold on
        %geoshow('landareas.shp', 'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'black');
        %axis off

        %SC1=scatterm(lat1(cruise==cruise_name(1)),lon1(cruise==cruise_name(1)),15,'s','filled','markerfacecolor',cmap(1,:),'markeredgecolor','none');
        %SC2=scatterm(lat1(cruise==cruise_name(2)),lon1(cruise==cruise_name(2)),15,'s','filled','markerfacecolor',cmap(2,:),'markeredgecolor','none');
        %SC3=scatterm(lat1(cruise==cruise_name(3)),lon1(cruise==cruise_name(3)),15,'s','filled','markerfacecolor',cmap(3,:),'markeredgecolor','none');
        %SC4=scatterm(lat1(cruise==cruise_name(4)),lon1(cruise==cruise_name(4)),15,'s','filled','markerfacecolor',cmap(4,:),'markeredgecolor','none');
        %SC5=scatterm(lat1(cruise==cruise_name(5)),lon1(cruise==cruise_name(5)),15,'s','filled','markerfacecolor',cmap(10,:),'markeredgecolor','none');
        %SC6=scatterm(lat1(cruise==cruise_name(6)),lon1(cruise==cruise_name(6)),15,'s','filled','markerfacecolor',cmap(12,:),'markeredgecolor','none');
        %SC7=scatterm(lat1(cruise==cruise_name(7)),lon1(cruise==cruise_name(7)),15,'s','filled','markerfacecolor',cmap(7,:),'markeredgecolor','none');

        %leg=legend([SC1,SC2,SC3,SC4,SC5,SC6,SC7],...
        %    {'Cruise 1','Cruise 2','Cruise 3','Cruise 4','Cruise 5','Cruise 6','Cruise 7'},'fontsize',8);


        %leg.ItemTokenSize(1) = 10;
        %leg.Location='southoutside';
        %leg.NumColumns=4;
        % sgtitle('POPCORN transects','FontSize', 20)

        %%

        str=strings([1,length(cruise_name)]);
        for i=1:length(cruise_name)
            str(i)=append('cruise',num2str(i));
        end

        % cmap(month(j),:); where j are the data points
        monthmap=  [0, 93, 174; 59, 63, 126; 169, 223, 205; 250, 250, 110;... % jan,feb same color
            252, 150, 176; 199, 61, 92; 177, 48, 104; 135, 73, 163;...
            204, 183, 255; 149, 167, 255; ...
            122, 186, 242; 0, 116, 217]./255;

        lat=myData_final(:,3);
        lon=myData_final(:,4);
        latlim=[min(lat),max(lat)];
        lonlim=[min(lon),max(lon)];
        POCavg_uM=myData_final(:,6);
        Pbiom_model_final=myData_final(:,9);
        month=myData_final(:,5);

        % FIGURE : NUMvsPOPCORN
        nexttile

        set(gcf,'color','w');

        scatter(POCavg_uM,Pbiom_model_final,[],'markerfacecolor','k');
        hold on
        t = linspace(0, 200, 10);
        plot(t,t,'k-','LineWidth',2)
        xlim([5 250])
        ylim([5 250])
        %xlabel('POC_{POPCORN} (\mugC/l)')
        %ylabel('POC_{NUM} (\mugC/l)')
        title('POC')
        set(gca,'xscale','log','yscale','log')

        ix = ~isnan(Pbiom_model_final);
        fit = (Pbiom_model_final(ix) ./ (POCavg_uM(ix)));
        fit = mean( fit(~isinf(fit)));

        %%
        % DETAILED Figure with latitudes



        nexttile

        for j=cruise1_ix
            %set(gca,'XTickLabel',[]);
            hold on
            plot(POCavg_uM(j),lat(j),'d','markerfacecolor','c','markeredgecolor','c','markersize',1);
            plot(Pbiom_model_final(j),lat(j),'d','markerfacecolor','k','markeredgecolor','k','markersize',1)
            xlim([0 250])
            ylim(latlim)
            %ylabel('latitude')
            %title(str(1))
        end

        %nexttile
        for j=cruise2_ix
            %set(gca,'XTickLabel',[]);
            %set(gca,'YTickLabel',[]);
            hold on
            plot(POCavg_uM(j),lat(j),'d','markerfacecolor','c','markeredgecolor','c','markersize',1);
            plot(Pbiom_model_final(j),lat(j),'d','markerfacecolor','k','markeredgecolor','k','markersize',1);
            xlim([0 250])
            ylim(latlim)
            %title(str(2))
        end

        %nexttile
        for j=cruise3_ix
            %set(gca,'XTickLabel',[]);
            %set(gca,'YTickLabel',[]);
            hold on
            plot(POCavg_uM(j),lat(j),'d','markerfacecolor','c','markeredgecolor','c','markersize',1);
            plot(Pbiom_model_final(j),lat(j),'d','markerfacecolor','k','markeredgecolor','k','markersize',1);
            xlim([0 250])
            ylim(latlim)
            %title(str(3))
        end

        %nexttile
        for j=cruise4_ix
            %set(gca,'XTickLabel',[]);
            %set(gca,'YTickLabel',[]);

            hold on
            plot(POCavg_uM(j),lat(j),'d','markerfacecolor','c','markeredgecolor','c','markersize',1);
            plot(Pbiom_model_final(j),lat(j),'d','markerfacecolor','k','markeredgecolor','k','markersize',1);
            xlim([0 250])
            ylim(latlim)
            %title(str(4))
        end


        %nexttile
        cruise5_ix_mooutlire=cruise5_ix;
        cruise5_ix_mooutlire(467)=[];
        for j=cruise5_ix_mooutlire
            % set(gca,'XTickLabel',[]);
            % set(gca,'YTickLabel',[]);
            hold on
            plot(POCavg_uM(j),lat(j),'d','markerfacecolor','c','markeredgecolor','c','markersize',1);
            plot(Pbiom_model_final(j),lat(j),'d','markerfacecolor','k','markeredgecolor','k','markersize',1);
            xlim([0 250])
            ylim(latlim)
            xlabel("biomass (\mugC/L)")
            %ylabel('latitude')
            %title(str(5))
        end

        %nexttile
        for j=cruise6_ix
            %set(gca,'YTickLabel',[]);
            hold on
            plot(POCavg_uM(j),lat(j),'d','markerfacecolor','c','markeredgecolor','c','markersize',1);
            plot(Pbiom_model_final(j),lat(j),'d','markerfacecolor','k','markeredgecolor','k','markersize',1);
            xlim([0 250])
            ylim(latlim)
            xlabel("biomass (\mugC/l)")
            %title(str(6))
        end

        %nexttile
        for j=cruise7_ix
            %set(gca,'YTickLabel',[]);
            hold on
            plot(POCavg_uM(j),lat(j),'d','markerfacecolor','c','markeredgecolor','c','markersize',1);
            plot(Pbiom_model_final(j),lat(j),'d','markerfacecolor','k','markeredgecolor','k','markersize',1);
            xlim([0 250])
            ylim(latlim)
            % ylabel('latitude')
            xlabel("biomass (\mugC/L)")
            %title(str(7))
        end

        %nexttile
        for j=cruise7_ix
            % set(gca,'YTickLabel',[]);
            hold on
            plot(POCavg_uM(j),lat(j),'d','markerfacecolor','c','markeredgecolor','c','markersize',1);
            plot(Pbiom_model_final(j),lat(j),'d','markerfacecolor','k','markeredgecolor','k','markersize',1);
            xlim([0 250])
            % ylim(latlim)
            xlabel("biomass (\mugC/L)")
            %title(str(7))
        end

        ylim(ylimLatitude)
    end
    function fit = PicophytoFiguresAtlantic(sim,r)
        load('poco_cpf_db_v2.mat','pococpfdbv2')
        load('pico_insitudata.mat','insitudata')

        %%
        % convert data variables to vectors
        datapoco=pococpfdbv2;
        insitu=table2array(datapoco(:,5));
        experiment=table2array(datapoco(:,15));

        % for i=1:length(insitudata)
        %     if (insitudata(i,2)>=30.63) && (insitudata(i,2)<=43.4)...
        %         && (insitudata(i,3)>=-10.33) && (insitudata(i,3)<=21.77)
        %        insitudata(i,:)=[];
        %     end
        % end


        % Assuming insitudata is a matrix with at least 3 columns (rows x columns)

        % Define latitude and longitude limits
        latLimits = [30.63, 43.4];
        lonLimits = [-10.33, 21.77];
        lonAtlantic = [-60, 22];
        % Create logical indices for rows that meet the criteria
        indicesToKeep = (insitudata(:, 2) >= latLimits(1) & insitudata(:, 2) <= latLimits(2)) & ...
            (insitudata(:, 3) >= lonLimits(1) & insitudata(:, 3) <= lonLimits(2));

        indicesAtlantic = (insitudata(:, 3) >= lonAtlantic(1) & insitudata(:, 3) <= lonAtlantic(2));

        % Use logical indexing to keep only the rows that meet the criteria
        filteredData = insitudata((indicesToKeep==false & indicesAtlantic), :);
        %%
        insitudata=filteredData;
        insitu=insitudata(:,1);
        lat=insitudata(:,2);
        lon=insitudata(:,3);
        month_poco=insitudata(:,6);

        % MODEL OUTPUT
        lonWrapped = wrapTo180(sim.x); % convert longitude ordinates from [0,360]->[-180,180]
        % r=rGroups;
        r_pico=find(r<=2);%-sim.p.nNutrients;
        r_pico=r_pico(4:end);
        Bpico1=sim.B(:,:,:,:,r_pico); % mugC/L
        Bpico_sum=squeeze(sum(Bpico1,5)); % sum all size-classes

        Bpico2=sim.Bpico;
        Bpico=Bpico_sum;
        %chlVolume_model=sim.ChlVolume(:,:,:,1); % top layer
        %%
        %------------------------------------------------------------------------------
        % Find model Biomass and Chla for the same coordinates and months as the data
        %------------------------------------------------------------------------------

        Picobiom_model=zeros(length(lat),1);
        chla_model=zeros(length(lat),1);

        for i=1:length(lat)

            [ ~, idx_lonG] = min( abs( lonWrapped-lon(i) ) );

            [ ~, idx_latG] = min( abs( sim.y-lat(i) ) );

            ix_lon(i)=idx_lonG;
            ix_lat(i)=idx_latG;

            % next step is to find combined coordinates
            Picobiom_model(i)=Bpico(month_poco(i),idx_lonG(1),idx_latG(1),1); % values taken at top layer
            %    chla_model(i,:)=chlVolume_model(month_poco(i),idx_lonG(1),idx_latG(1));
        end




        % FIGURE : NUM vs insitu
        nexttile
        set(gcf,'color','w');

        scatter(insitu,Picobiom_model,[],'markerfacecolor','k');
        % xlim([0 250])
        xlabel('')
        ylabel('Model')
        title('Pico')
        hold on
        plot([1 100],[1 100],'k-','linewidth',2)
        set(gca,'xscale','log','yscale','log')

        ix = ~isnan(Picobiom_model);
        fit = (Picobiom_model(ix) ./ insitu(ix));
        fit = mean( fit(~isinf(fit)));

        nexttile

        for j=1:length(lat)
            if (lon(j)>=-60 && lon(j)<=0)
                % set(gca,'XTickLabel',[]);
                hold on
                scatter(insitu(j),lat(j),[],'cyan','filled');
                scatter(Picobiom_model(j),lat(j),'kd','filled')
                % ylim(latlim)
                xlabel('Biomass (mg/m^3)')
                ylabel('Latitude')
                % title(str(1))
            end
        end
        legend('in situ','NUM')
        %title('Picophytoplankton biomass: insitu vs NUM')

        ylim(ylimLatitude)

    end

%
% Compare modelled copepod data along the AMT track with data from López and Anadón (2008)
% (As in Serra-Pompei et al 2022 figure 4).
%
    function [fit, Bmodel] = AMTcopepods(sim)
        %
        % Define data:
        %
        lat = [48.3178   47.2274   43.0218   40.0623   38.1931   34.7664   30.8723   26.1994   21.9938   20.7477   18.0997  10.0000    6.2617    2.0561   -6.5109  -10.4050  -14.6106  -18.8162  -22.5545  -26.4486  -29.5639  -32.679 -38.2866  -40.9346];
        lon = [-9.8450  -15.1765  -19.6695  -19.9979  -24.8281  -23.1558  -20.9826  -20.9748  -20.8011  -18.2991  -18.4613 -21.7812  -22.9416  -24.4346  -24.9203  -24.9138  -24.9068  -24.8998  -24.8936  -24.8871  -27.2152  -30.8767 -38.2007  -41.6963];

        % gC/m2:
        B(:,4) = [0.3773         0    0.5633    0.2040    0.2184    0.4332    0.1227    0.2330    0.2388    0.7008    0.7034    0.4824    0.4041    0.8479    0.2901    0.1400    0.0940   0.0608    0.1280    0.1464    0.0971    0.0146    2.0710    1.0978];
        B(:,3) = [0.3867         0    1.3984    0.1426    0.4378    0.2779    0.1462    0.1762    0.2137    0.1188    0.9154    0.3385    0.1510    0.4508    0.2764    0.0832    0.1096   0.0713    0.0841    0.0914    0.1161    0.0185    0.9881    2.2958];
        B(:,2) = [0.4546         0    3.3248    0.3611    0.6092    0.4867    0.4421    0.4759    0.7146    0.3288    1.0946    0.4472    0.3366    0.7643    0.2680    0.1965    0.2758   0.2222    0.2398    0.2159    0.2563    0.2960    0.9183    2.0695];
        B(:,1) = 0.001*[31.8556        0  334.3384  179.7190  474.8743  178.9145  171.5795   94.4342   76.2060  156.2178  323.9796  192.1518  171.0234   82.9817   52.5360   87.6493  134.3682 108.6737   91.5219   78.9438  127.7676  101.3492  234.6723  567.7467];
        Btot = sum(B,2);
        %
        % Extract copepods from the model:
        %
        for i = 1:length(lat)
            try
                ixWC = calcGlobalWatercolumn(lat(i), lon(i), sim);

                Bmodel(i) = 0;
                for j = find(sim.p.typeGroups>=10 & sim.p.typeGroups<=11)
                    ixB = sim.p.ixStart(j):sim.p.ixEnd(j) - sim.p.idxB+1;
                    B = sum(sim.B(:, ixWC.x, ixWC.y, :, ixB),5) / 1000;
                    B = squeeze(mean(B(sim.t>sim.t(end)-365,:,:,:),1));
                    ixZ = ~isnan(B);
                    B = sum(B(ixZ).*sim.dznom(ixZ));
                    Bmodel(i) = Bmodel(i) + B;
                end
            catch
                Bmodel(i) = NaN;
            end
        end
        %
        % Plot
        %
        nexttile
        loglog(Btot, Bmodel,'d','markerfacecolor','k','markeredgecolor','b')
        hold on
        plot([0.05 50],[0.05 50],'k-','linewidth',2)
        hold off
        mx = max(max(Btot), max(Bmodel));
        xlim([0 mx])
        ylim([0 mx])
        %xlabel('Data (g_C/m^2)')
        %ylabel('Model (g_C/m^2)')
        title('Copepods')

        ix = ~isnan(Bmodel);
        fit = (Bmodel(ix) ./ Btot(ix));
        fit = mean( fit(~isinf(fit)));

        nexttile
        plot(Btot,lat,'d','markerfacecolor','c','markeredgecolor','c')
        hold on
        plot(Bmodel, lat,'d','markerfacecolor','k','markeredgecolor','b')
        %ylabel('Latitude')
        ylim(ylimLatitude)
        xlabel('Copepod biomass (g_C/m^2)')
        %legend({'Data','Model'})


    end

    function NPPcafeVsNUMfinal(sim)
        load('npp_annual_avg_correct.mat','npp_annual_avg')
        load('NPPavg','nppAvg')

        lat=zeros(1080,1);
        lat(1)=1/6;
        for i=2:1080
            lat(i)=lat(i-1)+1/6;
        end
        lat=lat-90;
        %
        lon=zeros(2160,1);
        lon(1)=1/6;
        for i=2:2160
            lon(i)=lon(i-1)+1/6;
        end

        % nppAvgAnnual=squeeze(mean(nppAvg(1:11,:,:),1));
        nppCAFEinterpRec=interpolateToNUM(sim, npp_annual_avg);

        %%





        nexttile
        % cbar = panelGlobal(sim.x,sim.y,(sim.ProdNetAnnual(end,:,:)),[0,1000],...
        cbar = panelGlobal(sim.x,sim.y,(sim.ProdNetAnnual),[0,1000],...
            sTitle='Net primary production NUM', sProjection='eckert4');
        cbar.Label.String = '(mgC m^{-2}d^{-1})';
        cbar.Visible='on';
        % caxis([1,3])
        colorbar
        cmocean('haline',500)
        set(gca,'YTickLabel',[]);
        set(gca,'XTickLabel',[]);

        nexttile
        cbar = panelGlobal(sim.x,sim.y,(nppCAFEinterpRec'),[0,1000],...
            sTitle='Net primary production CAFE', sProjection='eckert4');
        cbar.Label.String = '(mgC m^{-2}d^{-1})';
        cbar.Visible='on';
        % caxis([1,3])
        colorbar
        cmocean('haline',500)
        set(gca,'YTickLabel',[]);
        set(gca,'XTickLabel',[]);


        % create meshgrid with CAFE's coordinates
        [lon1,lat1]=meshgrid(lon,lat);

        %create meshgrid with NUM's coordinates
        [lon2,lat2]=meshgrid(sim.x,sim.y);

        % interpolate NPP CAFE dimensions to NUM's
        nppCAFEinterp=interp2(lon1,lat1,squeeze(npp_annual_avg),lon2,lat2);

        orig0=(sim.ProdNetAnnual(end,:,:)); %NUM NPP
        orig=squeeze(orig0)'; % flip NUM npp data

        nppCAFEinterpRecenter= [nppCAFEinterp(:,65:end),nppCAFEinterp(:,1:64)];
        % make land white
        for i=1:64
            for j=1:128
                if (nppCAFEinterpRecenter(i,j)==-9999)
                    nppCAFEinterpRecenter(i,j)=NaN;
                    orig(i,j)=NaN;
                end
            end
        end

    end


end