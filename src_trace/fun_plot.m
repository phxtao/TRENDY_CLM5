function fun_plot(var_name, plot_unit, data_path, time_series, year_start, year_end, lat_grid, var_globe_time_series_diff, var_lat_time_series_diff)
%% plot lat time figures

fig = plot(time_series, var_globe_time_series_diff, 'LineStyle', '-', 'Color', 'blue', 'LineWidth', 2);
set(gca, 'TickDir', 'out');
set(gca, 'TickLength', [0.005, 0.005]);
set(gca, 'TickLength', [0.005, 0.005]);
set(gca, 'FontSize', 15);
% set(gca, 'ytick', (-90:15:90))
% set(gca, 'ytickLabel', (-90:15:90));
% set(gca, 'xtick', (year_start:50:year_end))
% set(gca, 'xTickLabel', (year_start:50:year_end));

hold on
plot([year_start, year_end], [0, 0], 'LineStyle', '--', 'Color', 'black', 'LineWidth', 2);

xlim([year_start year_end]);
title('', 'FontSize', 20)
ylabel(['Changes of ', var_name, ' (', plot_unit, ')'], 'FontSize', 20)
xlabel('Year', 'FontSize', 20)
set(gca, 'xcolor', 'black', 'ycolor', 'black')

fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 7 7];

print([data_path, 'globe_time_series_diff_', var_name, '.tif'], '-dtiffn', '-r300')
close


%% plot lat time figures
% generate colormap
white2red = [ones(101, 1), (1:-0.01:0)', (1:-0.01:0)'];
blue2white = [(0:0.01:1)', (0:0.01:1)', ones(101, 1)];
map_color = [blue2white; white2red];

% load coastlines

% lat time changes
fig = imagesc('XData', time_series, 'YData', lat_grid, 'CData', var_lat_time_series_diff);
set(gca, 'TickDir', 'out');
set(gca, 'TickLength', [0.005, 0.005]);
set(gca, 'TickLength', [0.005, 0.005]);
set(gca, 'FontSize', 15);
set(gca, 'ytick', (-90:15:90))
set(gca, 'ytickLabel', (-90:15:90));
set(gca, 'xtick', (year_start:50:year_end))
set(gca, 'xTickLabel', (year_start:50:year_end));

set(fig,'AlphaData',~isnan(var_lat_time_series_diff));
h = colorbar('westoutside');
ylabel(h, ['Changes of ', var_name, ' (', plot_unit, ')'], 'FontSize', 20)

caxis([-1 1]*max(max(var_lat_time_series_diff, [], 'omitnan'), [], 'omitnan'))
colormap(map_color)
% hold on
% plot(coastlon, coastlat, 'k', 'LineWidth', 0.5);
xlim([year_start year_end]);
ylim([-90 90]);
title('', 'FontSize', 20)
set(gca, 'xcolor', 'black', 'ycolor', 'black')

fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 12 5];

print([data_path, 'lat_time_series_diff_', var_name, '.tif'], '-dtiffn', '-r300')
close

end