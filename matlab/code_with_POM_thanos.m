
%% LAtitude vs Time --> Total biomass in upper 200m 
mAdult = logspace(log10(0.2), log10(1000), 7);
mCount = 1;
%
% lat = [50, 30, 20 0 -20 -40];
% lon = [-10, -25, -20 -25 -25 -45];
% Atlantic
% lat = linspace(60, -60, 29);
% lon = ones(1,29) * (-25);
% Pacific 
lat = 58;% linspace(58, -58, 21);
lon = -170;%ones(1,21) * (-170);
BsumTOTAL_area2 = zeros(length(lat), length(365*2));
BsumTOTAL_area_GENELASTIS = zeros(length(lat), length(365*2));
for idM = 1:length(lat)
    BsumTot = 0;
    % With POM
    n = 10;
    nCopepods = 10;
    nPOM = 10;
    p = setupNUMmodel(mAdult, n,nCopepods,nPOM);
    p = parametersWatercolumn(p);
    p.tEnd = 365*10;
    setHTL(0.15, 1, true, true);
    sim = simulateWatercolumn(p, lat(idM),lon(idM));
    
zcellbottom = sim.z + 0.5*sim.dznom;
    [~, iDepth200] = min(abs(zcellbottom-200));
    BB0 = sim.B(1:iDepth200,:,:); % BB(8x80x1095)
    newTime = find(sim.t > sim.t(end) - 365*2);
    for idTime = 1:length(newTime)
        BB = BB0(:,:, newTime(idTime)); % BB(8x90x1) one year
        BB2 = ones(iDepth200, 90);
        for ixDepth = 1:iDepth200
            BB2(ixDepth,:) = BB(ixDepth,:) * sim.dznom(ixDepth); % BB2(8x90x1)
        end
        BB2 = squeeze(sum(BB2(1:iDepth200,:),1)); % sum of upper 200 m, full days BB2(90x1)
        mCount
        mCount = mCount + 1;
        % Generalists
        BsumTOTAL_area_GENELASTIS(idM,idTime) = sum(BB2(1:10)); %SOS!!!!! POM!!!!!!!
        % Copepods
        BsumTOTAL_area2(idM,idTime) = sum(BB2(11:80)); %SOS!!!!! POM!!!!!!!
    end
end
%%
figure
for i = 1:length(lat)
    % plot(Bc(i,:) , lat_new(i) * ones(1,length(Bc(i,:))), '--o', 'LineWidth', 1.2)
    y = lat(i) * ones(1,length(BsumTOTAL_area2(i,:)));
    % k = BsumTOTAL_area2(i,:);
    k = BsumTOTAL_area_GENELASTIS(i,:);
    % k(k <= 0) = 0.0001;
    surf([newTime(:) newTime(:)], [y(:) y(:)], [k(:) k(:)], ...  % Reshape and replicate data
        'FaceColor', 'none', ...    % Don't bother filling faces with color
        'EdgeColor', 'interp', ...  % Use interpolated color for edges
        'LineWidth', 10);            % Make a thicker line
    view(2)
    set(gca, 'FontSize', 12)
    h = colorbar;
    %  set(gca,'ColorScale','log')
    % caxis([0.01 105])
    h.Label.String = "Biomass (mgC m^{-2})";
    h.FontSize = 12;
    % caxis([1.56 28360])
    caxis([20 2.1741e+04])
    xlim([newTime(1), newTime(end)])
    % ylim([min(lat_new) max(lat_new)])
    % title(['Annual mean biomass (ugC L-1) of the upper 200 m, mHTL = ', num2str(mHTL_new), ' ugC, mortHTL = ', num2str(mortHTL_new) ])
    ylabel('Latitude (degrees)')
    datetick('x', 'mmm', 'keeplimits')
    xlabel('Time (days)')
    set(gca,'ColorScale','log')
    hold on
end
% title('Generalists')

