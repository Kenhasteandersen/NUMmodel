                
figure(1)
clf(1)
set(gcf,'color','w');

for i = 1:length(newRemin2)

    for j = 1:length(newSinkingPOM)

        plot(NPP_cell_month_mean{i, j}(end-11:end), 'o-r', 'LineWidth', 1)% mgC/m2/day
        % plot(NPP_cell_month_mean{i, j}, 'o-r', 'LineWidth', 1)% mgC/m2/day
        hold on
        plot(NPP_extracted(1,:), 'ob--')
        plot(NPP_extracted(2,:), 'og--')
        plot(NPP_extracted(3,:), 'om--')

xlabel('Time (month)')
ylabel('NPP (mgC / m^2 /day)')
 mTitle = append('Lat: ', string(lat_to_find), ', Lon: ', string(lon_to_find));
 my_title = append('remin2 = ', string(newRemin2(i)), ' / sinking = ', string(newSinkingPOM(j)));
        title(my_title)

    end
end
title(t,mTitle, 'FontSize', 24)
lg  = legend('NUM','Eppley model', 'Standard VGPM', 'CAFE','NumColumns',2); 
lg.Layout.Tile = 'South'; % <-- Legend placement with tiled layout

%%
%  Compare NPP over the last years in the model        
set(gcf,'color','w');
figure(2)
clf(2)
        % total_months=size(matr_NPP,2;)  % last year= 4th year
        plot(NPP_cell_month_mean{i, j}(end-11:end), 'o-r', 'LineWidth', 1)% mgC/m2/day
        hold on                                     % 3rd year from 15years+
        plot(NPP_cell_month_mean{i, j}( (noYears-2)*12-1:(noYears-1)*12), 'o-g', 'LineWidth', 1)% mgC/m2/day
                                                     % 2nd year                       
        plot(NPP_cell_month_mean{i, j}( (noYears-3)*12-1:(noYears-2)*12), 'o-m', 'LineWidth', 1)% mgC/m2/day
        plot(NPP_cell_month_mean{i, j}( (noYears-4)*12-1:(noYears-3)*12), 'o-b', 'LineWidth', 1)% mgC/m2/day
legend('year 20','year 19','year 18','year 17')
xlabel('Time (month)')
ylabel('NPP (mgC / m^2 /day)')
axis tight







% Calculate phyoplankton 60% photosynthetic
 % Bph_orig=calcPhyto(sim,simR);
Bph=Bph_orig;
z = sim.z + 0.5*sim.dznom;
% Make a layer at z = 0 with the same value as in the first grid point:
%
t = sim.t;
z = [0; z];
Bph(:,2:length(z),:) = Bph;
Bph(:,1,:) = Bph(:,2,:);

% ylimit = [-max(z), 0];
ylimit = [-400, 0];
xlimit = [sim.t(1) sim.t(end)];
%
figure(3)
clf(3)
set(gcf,'Color','w');
tiledlayout(3,1)
nexttile
    % Bph(Bph < 0.01) = 0.01; % Set low biomasses to the lower limit to avoid white space in plot
    contourf(t,-z,log10(squeeze(sum(Bph(:,:,:),3)))',[linspace(-2,3,5)],'LineStyle','none')

    title('Phytoplankton _{0.6}')
    ylabel('Depth (m)')
    axis tight
    colorbar
    % caxis([0.1 100])
    colorbar('ticks',-2:3,'limits',[-2 3])
    % clim([0, 1000])
    ylim(ylimit)
    xlim(xlimit)
    xlabel('Time (days)')
    a = colorbar;
    a.Label.String = 'Concentartion log_{10}\mugC/l';

    nexttile
    % Bph(Bph < 0.01) = 0.01; % Set low biomasses to the lower limit to avoid white space in plot
    contourf(t(end-364:end),-z,log10(squeeze(sum(Bph(end-364:end,:,:),3)))',[linspace(-2,3,5)],'LineStyle','none')

    title('Phytoplankton _{0.6}  (Final year)')
    ylabel('Depth (m)')
    axis tight
    colorbar
    % caxis([0.1 100])
    colorbar('ticks',-2:3,'limits',[-2 3])
    % clim([0, 1000])
    ylim(ylimit)
    % xlim(xlimit)
    xlabel('Time (days)')
    a = colorbar;
    a.Label.String = 'Concentartion log_{10}\mugC/l';

    nexttile
    % Bph(Bph < 0.01) = 0.01; % Set low biomasses to the lower limit to avoid white space in plot
    contourf(t(end-2*365+1:end-365),-z,log10(squeeze(sum(Bph(end-2*365+1:end-365,:,:),3)))',[linspace(-2,3,5)],'LineStyle','none')

    title('Phytoplankton _{0.6}  (2nd to last year)')
    ylabel('Depth (m)')
    axis tight
    colorbar
    % caxis([0.1 100])
    colorbar('ticks',-2:3,'limits',[-2 3])
    % clim([0, 1000])
    ylim(ylimit)
    % xlim(xlimit)
    xlabel('Time (days)')
    a = colorbar;
    a.Label.String = 'Concentartion log_{10}\mugC/l';
