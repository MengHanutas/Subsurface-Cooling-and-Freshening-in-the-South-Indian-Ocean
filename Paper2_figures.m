%% Figure plot for paper 2, Figure 1
%% load data
load Argo_new.mat
Argo.dh=nanmean(Argo.geo_strf_dyn_height,4);
Argo.gm=nanmean(Argo.gamma_n,4);
load coastlines.mat


%calculate alpha and beta:
% Initialize matrices to store alpha and beta values
[nlon, nlat, ndepth, ntime] = size(Argo.CT_month_all);
alpha = nan(nlon, nlat, ndepth, ntime);
beta = nan(nlon, nlat, ndepth, ntime);

% Loop through time and depth
for t = 1:ntime
    for d = 1:ndepth
        % Extract sea pressure at this depth level (constant across lat, lon)
        p_level = Argo.pressure(d); % sea pressure for the specific depth in decibars
        
        % Compute alpha (thermal expansion coefficient)
        alpha(:, :, d, t) = gsw_alpha(Argo.SA_month_all(:, :, d, t), Argo.CT_month_all(:, :, d, t), p_level);
        
        % Compute beta (saline contraction coefficient)
        beta(:, :, d, t) = gsw_beta(Argo.SA_month_all(:, :, d, t), Argo.CT_month_all(:, :, d, t), p_level);
    end
end

%calculate alpha and beta mean
alpha_mean=nanmean(alpha,4);
beta_mean=nanmean(beta,4);

%% Paper 2 Figure 1

%%
%% ========== 1. Figure setup ==========
figure('Units', 'inches', 'Position', [1 1 8.5 7], 'PaperPositionMode', 'auto');

% Define margins and gaps (normalized units)
left_margin   = 0.06;
right_margin  = 0.08;   % Keep small; colorbars are inside subplots
bottom_margin = 0.09;
top_margin    = 0.07;   % Room for labels above subplots

% REDUCED GAPS – much tighter spacing
h_gap         = 0.03;  % horizontal gap between columns
v_gap         = 0.025; % vertical gap between rows

nrows = 2;
ncols = 3;

% Compute subplot width and height
sub_w = (1 - left_margin - right_margin - (ncols-1)*h_gap) / ncols;
sub_h = (1 - bottom_margin - top_margin - (nrows-1)*v_gap) / nrows;

% Pre-compute subplot positions: [x, y, width, height]
positions = zeros(nrows*ncols, 4);
for row = 1:nrows
    for col = 1:ncols
        idx = (row-1)*ncols + col;
        x = left_margin + (col-1)*(sub_w + h_gap);
        y = bottom_margin + (nrows-row)*(sub_h + v_gap);  % top row first
        positions(idx, :) = [x, y, sub_w, sub_h];
    end
end

% Load land and coastline data once (assumes Mapping Toolbox)
land = shaperead('landareas.shp', 'UseGeoCoords', true);
load coastlines

% ========== 2. Define depth indices and variables ==========
depth_ranges = {1:25, 25:44, 44:size(Argo.CT_month_all,3)};  % 0-300, 300-1000, 1000-2000 m
depth_labels = {'0-300 m', '300-1000 m', '1000-2000 m'};
var_names = {'CT', 'SA'};
var_units = {'[^{\circ}C/Yr]', '[g/kg/Yr]'};
caxis_limits = {[-0.04 0.04], [-0.005 0.005]};

% ========== 3. Loop over variables (CT and SA) and depth layers ==========
for var_idx = 1:2
    var_name = var_names{var_idx};
    for depth_idx = 1:3
        % Subplot index (row-major order)
        sub_idx = (var_idx-1)*3 + depth_idx;
        subplot('Position', positions(sub_idx, :));
        hold on;
        
        % 3a. Extract and deseason the data
        depth_inds = depth_ranges{depth_idx};
        data_3d = squeeze(mean(Argo.([var_name '_month_all'])(:, :, depth_inds, :), 3));
        data_deseason = deseason(data_3d, Argo.t2, 'monthly');
        
        % 3b. Compute linear trend (per year) and significance (ESS corrected)
        [trend_map, p_map] = calc_trend_ess(data_deseason);
        sig_mask = p_map >= 0.05;   % non-significant areas
        
        % 3c. Plot trend as pcolor
        pcolor(Argo.longitude, Argo.latitude, trend_map');
        shading flat;
        
        % 3d. Stipple insignificant areas
        stipple(Argo.longitude, Argo.latitude, sig_mask', ...
                'color', 'k', 'density', 50, 'marker', '.', 'markersize', 4);
        
        % 3e. Overlay dynamic height contours
        dh_mean = squeeze(mean(Argo.dh(:, depth_inds, :), 2));
        [C, h] = contour(Argo.longitude, Argo.latitude, dh_mean', 'k-', 'LineWidth', 0.5);
        clabel(C, h, 'FontSize', 7, 'Color', 'k');
        
        % 3f. Add land mask (grey filled) and coastline
        geoshow(land, 'FaceColor', [0.85 0.85 0.85], 'EdgeColor', 'none');
        plot(coastlon, coastlat, 'k-', 'LineWidth', 0.6);
        
        % 3g. Set colormap and color axis
        cmap_full = cmocean('balance', 36);
        cmap_full = cmap_full([1:16, 21:36], :);   % 15 blue + 15 red, perfectly balanced
        colormap(gca, cmap_full);
        caxis(caxis_limits{var_idx});
        
        % 3h. Add subplot label (a, b, c, ...) at top-left corner
        label_char = char(96 + sub_idx);  % a=97, b=98, ...
        text(0.02, 1.08, ['(' label_char ')'], 'Units', 'normalized', ...
             'FontSize', 11, 'FontWeight', 'bold', 'Color', 'k', ...
             'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');
        
        % 3i. Axis labels and formatting using apply_geo_ticks
        showX = (sub_idx >= 4);                 % Only bottom row gets X labels
        showY = (sub_idx == 1 || sub_idx == 4); % Only left column gets Y labels
        
        set(gca, 'tickdir', 'out', 'FontSize', 9);
        xlim([20 150]);
        ylim([-65 0]);
        
        % Apply custom geographic tick formatting
        apply_geo_ticks(gca, 20:20:140, -60:20:0, showX, showY);
        
        if showX
            xlabel('Longitude');
        end
        if showY
            ylabel('Latitude');
        end
        
        % 3j. Add colorbar only for the rightmost subplots (3rd column)
        if depth_idx == 3
            c = colorbar;
            c.Label.String = var_units{var_idx};
            c.TickDirection = 'out';  
            % If you want specific ticks (e.g., for salinity):
            % c.Ticks = -0.005:0.001:0.005; 
        end
        
        % 3k. Add depth label as text inside subplot
        text(0.98, 0.62, depth_labels{depth_idx}, 'Units', 'normalized', ...
             'FontSize', 8, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');
        
        hold off;
    end
end
%% Figure 2
%% ========== 1. Figure setup ==========
figure('Units', 'inches', 'Position', [1 1 8.5 7], 'PaperPositionMode', 'auto');

left_margin   = 0.06;
right_margin  = 0.08;
bottom_margin = 0.09;
top_margin    = 0.06;
h_gap         = 0.02;
v_gap         = 0.04;

nrows = 2;
ncols = 3;

sub_w = (1 - left_margin - right_margin - (ncols-1)*h_gap) / ncols;
sub_h = (1 - bottom_margin - top_margin - (nrows-1)*v_gap) / nrows;

positions = zeros(nrows*ncols, 4);
for row = 1:nrows
    for col = 1:ncols
        idx = (row-1)*ncols + col;
        x = left_margin + (col-1)*(sub_w + h_gap);
        y = bottom_margin + (nrows-row)*(sub_h + v_gap);
        positions(idx, :) = [x, y, sub_w, sub_h];
    end
end

% Load land and coastline data
land = shaperead('landareas.shp', 'UseGeoCoords', true);
load coastlines

% ---- Geographic tick settings (used by apply_geo_ticks) ----
lonTicks = 20:20:140;    % 20E, 40E, ..., 140E
latTicks = -60:20:0;     % 60S, 40S, 20S, 0

% ========== 2. Common plotting parameters ==========
caxis_lim = [-1e-5 1e-5];
cmap_full = cmocean('balance', 40);
cmap_modified = cmap_full([1:19, 22:40], :);

% ========== 3. Panel (a): alpha*T|z 0-300m ==========
ax = subplot('Position', positions(1,:));
t_z_ap = squeeze(nanmean(alpha_mean(:,:,1:25),3)) .* ...
         squeeze(nanmean(Argo.t_z(:,:,1:25,:),3));
t_z_ap_ds = deseason(t_z_ap, Argo.t2, 'monthly');
%[t_z_ap_tr, p] = trend(t_z_ap_ds, 12, 'omitnan');      
[t_z_ap_tr, p] = calc_trend_ess(t_z_ap_ds);
%sig_mask = p >= 0.05;

pcolor(Argo.longitude, Argo.latitude, t_z_ap_tr'); shading flat;
hold on;
stipple(Argo.longitude, Argo.latitude, (p>=0.05)', ...
        'color','k','density',50,'marker','.','markersize',5);
[c,h] = contour(Argo.longitude, Argo.latitude, ...
                squeeze(nanmean(Argo.gm(:,:,1:25),3))', ...
                'k-','LineWidth',0.5);
clabel(c,h,'FontSize',7,'Color','k');
geoshow(land,'FaceColor',[0.85 0.85 0.85],'EdgeColor','none');
plot(coastlon, coastlat, 'k-','LineWidth',0.6);

ylabel('Latitude');
set(ax,'tickdir','out','FontSize',9);
xlim([20 150]); ylim([-65 0]);
apply_geo_ticks(ax, lonTicks, latTicks, false, true);   % hide X, show Y

text(0.83, 0.6, '$\alpha T''|_{p}$','Units','normalized', ...
     'FontSize',11,'FontWeight','bold','Color','k','Interpreter','latex');

colormap(ax, cmap_modified); caxis(caxis_lim);
hold off;

% ========== 4. Panel (b): alpha*T spice 0-300m ==========
ax = subplot('Position', positions(2,:));
t_spice_ap = squeeze(nanmean(alpha_mean(:,:,1:25),3)) .* ...
             squeeze(nanmean(Argo.t_spice_z(:,:,1:25,:),3));
t_spice_ap_ds = deseason(fillmissing(t_spice_ap,'linear',3), Argo.t2,'monthly');
%[t_spice_ap_tr, p] = trend(t_spice_ap_ds, 12,'omitnan');
[t_spice_ap_tr, p] = calc_trend_ess(t_spice_ap_ds);

pcolor(Argo.longitude, Argo.latitude, t_spice_ap_tr'); shading flat;
hold on;
stipple(Argo.longitude, Argo.latitude, (p>=0.05)', ...
        'color','k','density',50,'marker','.','markersize',5);
[c,h] = contour(Argo.longitude, Argo.latitude, ...
                squeeze(nanmean(Argo.gm(:,:,1:25),3))', ...
                'k-','LineWidth',0.5);
clabel(c,h,'FontSize',7,'Color','k');
geoshow(land,'FaceColor',[0.85 0.85 0.85],'EdgeColor','none');
plot(coastlon, coastlat, 'k-','LineWidth',0.6);

set(ax,'tickdir','out','FontSize',9);
xlim([20 150]); ylim([-65 0]);
apply_geo_ticks(ax, lonTicks, latTicks, false, false);  % hide X and Y

text(0.83, 0.6, '$\alpha T''|_{\mathrm{spice}}$','Units','normalized', ...
     'FontSize',11,'FontWeight','bold','Color','k','Interpreter','latex');

colormap(ax, cmap_modified); caxis(caxis_lim);
hold off;

% ========== 5. Panel (c): alpha*T heave 0-300m ==========
ax = subplot('Position', positions(3,:));
t_heave_ap = squeeze(nanmean(alpha_mean(:,:,1:25),3)) .* ...
             squeeze(nanmean(Argo.t_heave_z(:,:,1:25,:),3));
t_heave_ap_ds = deseason(fillmissing(t_heave_ap,'linear',3), Argo.t2,'monthly');
%[t_heave_ap_tr, p] = trend(t_heave_ap_ds, 12,'omitnan');
[t_heave_ap_tr, p] = calc_trend_ess(t_heave_ap_ds);

pcolor(Argo.longitude, Argo.latitude, -t_heave_ap_tr'); shading flat;
hold on;
stipple(Argo.longitude, Argo.latitude, (p>=0.05)', ...
        'color','k','density',50,'marker','.','markersize',5);
[c,h] = contour(Argo.longitude, Argo.latitude, ...
                squeeze(nanmean(Argo.gm(:,:,1:25),3))', ...
                'k-','LineWidth',0.5);
clabel(c,h,'FontSize',7,'Color','k');
geoshow(land,'FaceColor',[0.85 0.85 0.85],'EdgeColor','none');
plot(coastlon, coastlat, 'k-','LineWidth',0.6);

c = colorbar; c.Label.String = '[kg/m^3/Yr]'; c.TickDirection = 'out';  
set(ax,'tickdir','out','FontSize',9);
xlim([20 150]); ylim([-65 0]);
apply_geo_ticks(ax, lonTicks, latTicks, false, false);  % hide X and Y

text(0.83, 0.6, '$\alpha T''|_{\mathrm{heave}}$','Units','normalized', ...
     'FontSize',11,'FontWeight','bold','Color','k','Interpreter','latex');

colormap(ax, cmap_modified); caxis(caxis_lim);
hold off;

% ========== 6. Panel (d): beta*S|z 0-300m ==========
ax = subplot('Position', positions(4,:));
s_z_bt = squeeze(nanmean(beta_mean(:,:,1:25),3)) .* ...
         squeeze(nanmean(Argo.s_z(:,:,1:25,:),3));
s_z_bt_ds = deseason(s_z_bt, Argo.t2,'monthly');
%[s_z_bt_tr, p] = trend(s_z_bt_ds, 12,'omitnan');
[s_z_bt_tr, p] = calc_trend_ess(s_z_bt_ds);

pcolor(Argo.longitude, Argo.latitude, s_z_bt_tr'); shading flat;
hold on;
stipple(Argo.longitude, Argo.latitude, (p>=0.05)', ...
        'color','k','density',50,'marker','.','markersize',5);
[c,h] = contour(Argo.longitude, Argo.latitude, ...
                squeeze(nanmean(Argo.gm(:,:,1:25),3))', ...
                'k-','LineWidth',0.5);
clabel(c,h,'FontSize',7,'Color','k');
geoshow(land,'FaceColor',[0.85 0.85 0.85],'EdgeColor','none');
plot(coastlon, coastlat, 'k-','LineWidth',0.6);

ylabel('Latitude');
xlabel('Longitude');
set(ax,'tickdir','out','FontSize',9);
xlim([20 150]); ylim([-65 0]);
apply_geo_ticks(ax, lonTicks, latTicks, true, true);    % show X and Y

text(0.83, 0.6, '$\beta S''|_{p}$','Units','normalized', ...
     'FontSize',11,'FontWeight','bold','Color','k','Interpreter','latex');

colormap(ax, cmap_modified); caxis(caxis_lim);
hold off;

% ========== 7. Panel (e): beta*S spice 0-300m ==========
ax = subplot('Position', positions(5,:));
s_spice_bt = squeeze(nanmean(beta_mean(:,:,1:25),3)) .* ...
             squeeze(nanmean(Argo.s_spice_z(:,:,1:25,:),3));
s_spice_bt_ds = deseason(fillmissing(s_spice_bt,'linear',3), Argo.t2,'monthly');
%[s_spice_bt_tr, p] = trend(s_spice_bt_ds, 12,'omitnan');
[s_spice_bt_tr, p] = calc_trend_ess(s_spice_bt_ds);

pcolor(Argo.longitude, Argo.latitude, s_spice_bt_tr'); shading flat;
hold on;
stipple(Argo.longitude, Argo.latitude, (p>=0.05)', ...
        'color','k','density',50,'marker','.','markersize',5);
[c,h] = contour(Argo.longitude, Argo.latitude, ...
                squeeze(nanmean(Argo.gm(:,:,1:25),3))', ...
                'k-','LineWidth',0.5);
clabel(c,h,'FontSize',7,'Color','k');
geoshow(land,'FaceColor',[0.85 0.85 0.85],'EdgeColor','none');
plot(coastlon, coastlat, 'k-','LineWidth',0.6);

xlabel('Longitude');
set(ax,'tickdir','out','FontSize',9);
xlim([20 150]); ylim([-65 0]);
apply_geo_ticks(ax, lonTicks, latTicks, true, false);   % show X, hide Y

text(0.83, 0.6, '$\beta S''|_{\mathrm{spice}}$','Units','normalized', ...
     'FontSize',11,'FontWeight','bold','Color','k','Interpreter','latex');

colormap(ax, cmap_modified); caxis(caxis_lim);
hold off;

% ========== 8. Panel (f): beta*S heave 0-300m ==========
ax = subplot('Position', positions(6,:));
s_heave_bt = squeeze(nanmean(beta_mean(:,:,1:25),3)) .* ...
             squeeze(nanmean(Argo.s_heave_z(:,:,1:25,:),3));
s_heave_bt_ds = deseason(fillmissing(s_heave_bt,'linear',3), Argo.t2,'monthly');
%[s_heave_bt_tr, p] = trend(-s_heave_bt_ds, 12,'omitnan');
[s_heave_bt_tr, p] = calc_trend_ess(-s_heave_bt_ds);

pcolor(Argo.longitude, Argo.latitude, s_heave_bt_tr'); shading flat;
hold on;
stipple(Argo.longitude, Argo.latitude, (p>=0.05)', ...
        'color','k','density',50,'marker','.','markersize',5);
[c,h] = contour(Argo.longitude, Argo.latitude, ...
                squeeze(nanmean(Argo.gm(:,:,1:25),3))', ...
                'k-','LineWidth',0.5);
clabel(c,h,'FontSize',7,'Color','k');
geoshow(land,'FaceColor',[0.85 0.85 0.85],'EdgeColor','none');
plot(coastlon, coastlat, 'k-','LineWidth',0.6);

c = colorbar; c.Label.String = '[kg/m^3/Yr]'; c.TickDirection = 'out';  
xlabel('Longitude');
set(ax,'tickdir','out','FontSize',9);
xlim([20 150]); ylim([-65 0]);
apply_geo_ticks(ax, lonTicks, latTicks, true, false);   % show X, hide Y

text(0.83, 0.6, '$\beta S''|_{\mathrm{heave}}$','Units','normalized', ...
     'FontSize',11,'FontWeight','bold','Color','k','Interpreter','latex');

colormap(ax, cmap_modified); caxis(caxis_lim);
hold off;

%========== 9. Add (a), (b), ... labels above each subplot ==========
for idx = 1:nrows*ncols
    x_left = positions(idx, 1);
    y_top  = positions(idx, 2) + positions(idx, 4);
    label_str = ['(' char(96+idx) ')'];
    annotation('textbox', [x_left + 0.01, y_top + 0.04, 0, 0], ...
               'String', label_str, ...
               'FitBoxToText', 'on', ...
               'EdgeColor', 'none', ...
               'FontSize', 12, ...
               'FontWeight', 'bold');
end

% ========== Local helper functions ==========
function apply_geo_ticks(ax, lonT, latT, showX, showY)
% APPLY_GEO_TICKS  Degree-formatted lon/lat tick labels.
    if nargin < 4 || isempty(showX), showX = true; end
    if nargin < 5 || isempty(showY), showY = true; end

    set(ax, 'XTick', lonT, 'YTick', latT);

    if showX
        set(ax, 'XTickLabel', ...
            arrayfun(@lon_label, lonT, 'UniformOutput', false));
    else
        set(ax, 'XTickLabel', []);
    end
    if showY
        set(ax, 'YTickLabel', ...
            arrayfun(@lat_label, latT, 'UniformOutput', false));
    else
        set(ax, 'YTickLabel', []);
    end
end

function lbl = lon_label(v)
    v = round(v);
    if v == 0
        lbl = '0^{\circ}';
    elseif v > 0
        lbl = sprintf('%d\\circE', v);
    else
        lbl = sprintf('%d\\circW', abs(v));
    end
end

function lbl = lat_label(v)
    v = round(v);
    if v == 0
        lbl = '0^{\circ}';
    elseif v > 0
        lbl = sprintf('%d\\circN', v);
    else
        lbl = sprintf('%d\\circS', abs(v));
    end
end
% ========== 10. Optional: Add depth label inside each panel ==========
% (Uncomment if you want to show depth range inside each subplot)
% for idx = 1:3
%     subplot('Position', positions(idx,:));
%     text(0.98, 0.02, '0-300 m', 'Units', 'normalized', ...
%          'FontSize', 8, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');
% end


%% Figure 3
%% ========== 1. Figure setup ==========
figure('Units', 'inches', 'Position', [1 1 8.5 7], 'PaperPositionMode', 'auto');

left_margin   = 0.06;
right_margin  = 0.08;
bottom_margin = 0.1;
top_margin    = 0.06;
h_gap         = 0.02;
v_gap         = 0.04;

nrows = 2;
ncols = 3;

sub_w = (1 - left_margin - right_margin - (ncols-1)*h_gap) / ncols;
sub_h = (1 - bottom_margin - top_margin - (nrows-1)*v_gap) / nrows;

positions = zeros(nrows*ncols, 4);
for row = 1:nrows
    for col = 1:ncols
        idx = (row-1)*ncols + col;
        x = left_margin + (col-1)*(sub_w + h_gap);
        y = bottom_margin + (nrows-row)*(sub_h + v_gap);
        positions(idx, :) = [x, y, sub_w, sub_h];
    end
end

land = shaperead('landareas.shp', 'UseGeoCoords', true);
load coastlines

% ---- Geographic tick settings ----
lonTicks = 20:20:140;    % 20E ... 140E
latTicks = -60:20:0;     % 60S ... 0

% ========== 2. Common plotting parameters ==========
caxis_lim = [-5e-6 5e-6];
cmap_full = cmocean('balance', 36);
cmap_modified = cmap_full([1:16, 21:36], :);

% ========== 3. Panel (a): alpha*T|z 300-1000m ==========
ax = subplot('Position', positions(1,:));
t_z_ap_300_1000 = squeeze(nanmean(alpha_mean(:,:,25:44),3)) .* ...
                  squeeze(nanmean(Argo.t_z(:,:,25:44,:),3));
t_z_ap_ds = deseason(t_z_ap_300_1000, Argo.t2, 'monthly');
%[t_z_ap_tr, p] = trend(t_z_ap_ds, 12, 'omitnan');
[t_z_ap_tr, p] = calc_trend_ess(t_z_ap_ds);

pcolor(Argo.longitude, Argo.latitude, t_z_ap_tr'); shading flat;
hold on;
stipple(Argo.longitude, Argo.latitude, (p>=0.05)', ...
        'color','k','density',50,'marker','.','markersize',5);
[c,h] = contour(Argo.longitude, Argo.latitude, ...
                squeeze(nanmean(Argo.gm(:,:,25:44),3))', ...
                'k-','LineWidth',0.5);
clabel(c,h,'FontSize',7,'Color','k');
geoshow(land,'FaceColor',[0.85 0.85 0.85],'EdgeColor','none');
plot(coastlon, coastlat, 'k-','LineWidth',0.6);

ylabel('Latitude');
set(ax,'tickdir','out','FontSize',9);
xlim([20 150]); ylim([-65 0]);
apply_geo_ticks(ax, lonTicks, latTicks, false, true);   % hide X, show Y

text(0.83, 0.6, '$\alpha T''|_{p}$','Units','normalized', ...
     'FontSize',11,'FontWeight','bold','Color','k','Interpreter','latex');

colormap(ax, cmap_modified); caxis(caxis_lim);
hold off;

% ========== 4. Panel (b): alpha*T spice 300-1000m ==========
ax = subplot('Position', positions(2,:));
t_spice_ap_300_1000 = squeeze(nanmean(alpha_mean(:,:,25:44),3)) .* ...
                      squeeze(nanmean(Argo.t_spice_z(:,:,25:44,:),3));
t_spice_ap_ds = deseason(fillmissing(t_spice_ap_300_1000,'linear',3), ...
                         Argo.t2, 'monthly');
%[t_spice_ap_tr, p] = trend(t_spice_ap_ds, 12, 'omitnan');
[t_spice_ap_tr, p] = calc_trend_ess(t_spice_ap_ds);

pcolor(Argo.longitude, Argo.latitude, t_spice_ap_tr'); shading flat;
hold on;
stipple(Argo.longitude, Argo.latitude, (p>=0.05)', ...
        'color','k','density',50,'marker','.','markersize',5);
[c,h] = contour(Argo.longitude, Argo.latitude, ...
                squeeze(nanmean(Argo.gm(:,:,25:44),3))', ...
                'k-','LineWidth',0.5);
clabel(c,h,'FontSize',7,'Color','k');
geoshow(land,'FaceColor',[0.85 0.85 0.85],'EdgeColor','none');
plot(coastlon, coastlat, 'k-','LineWidth',0.6);

set(ax,'tickdir','out','FontSize',9);
xlim([20 150]); ylim([-65 0]);
apply_geo_ticks(ax, lonTicks, latTicks, false, false);  % hide X and Y

text(0.83, 0.6, '$\alpha T''|_{\mathrm{spice}}$','Units','normalized', ...
     'FontSize',11,'FontWeight','bold','Color','k','Interpreter','latex');

colormap(ax, cmap_modified); caxis(caxis_lim);
hold off;

% ========== 5. Panel (c): alpha*T heave 300-1000m ==========
ax = subplot('Position', positions(3,:));
t_heave_ap_300_1000 = squeeze(nanmean(alpha_mean(:,:,25:44),3)) .* ...
                      squeeze(nanmean(Argo.t_heave_z(:,:,25:44,:),3));
t_heave_ap_ds = deseason(fillmissing(t_heave_ap_300_1000,'linear',3), ...
                         Argo.t2, 'monthly');
%[t_heave_ap_tr, p] = trend(t_heave_ap_ds, 12, 'omitnan');
[t_heave_ap_tr, p] = calc_trend_ess(t_heave_ap_ds);

pcolor(Argo.longitude, Argo.latitude, -t_heave_ap_tr'); shading flat;
hold on;
stipple(Argo.longitude, Argo.latitude, (p>=0.05)', ...
        'color','k','density',50,'marker','.','markersize',5);
[c,h] = contour(Argo.longitude, Argo.latitude, ...
                squeeze(nanmean(Argo.gm(:,:,25:44),3))', ...
                'k-','LineWidth',0.5);
clabel(c,h,'FontSize',7,'Color','k');
geoshow(land,'FaceColor',[0.85 0.85 0.85],'EdgeColor','none');
plot(coastlon, coastlat, 'k-','LineWidth',0.6);

c = colorbar; c.Label.String = '[kg/m^3/Yr]'; c.TickDirection = 'out'; c.Ticks = -5e-6:1e-6:5e-6;
set(ax,'tickdir','out','FontSize',9);
xlim([20 150]); ylim([-65 0]);
apply_geo_ticks(ax, lonTicks, latTicks, false, false);  % hide X and Y

text(0.83, 0.6, '$\alpha T''|_{\mathrm{heave}}$','Units','normalized', ...
     'FontSize',11,'FontWeight','bold','Color','k','Interpreter','latex');

colormap(ax, cmap_modified); caxis(caxis_lim);
hold off;

% ========== 6. Panel (d): beta*S|z 300-1000m ==========
ax = subplot('Position', positions(4,:));
s_z_bt_300_1000 = squeeze(nanmean(beta_mean(:,:,25:44),3)) .* ...
                  squeeze(nanmean(Argo.s_z(:,:,25:44,:),3));
s_z_bt_ds = deseason(s_z_bt_300_1000, Argo.t2, 'monthly');
%[s_z_bt_tr, p] = trend(s_z_bt_ds, 12, 'omitnan');
[s_z_bt_tr, p] = calc_trend_ess(s_z_bt_ds);

pcolor(Argo.longitude, Argo.latitude, s_z_bt_tr'); shading flat;
hold on;
stipple(Argo.longitude, Argo.latitude, (p>=0.05)', ...
        'color','k','density',50,'marker','.','markersize',5);
[c,h] = contour(Argo.longitude, Argo.latitude, ...
                squeeze(nanmean(Argo.gm(:,:,25:44),3))', ...
                'k-','LineWidth',0.5);
clabel(c,h,'FontSize',7,'Color','k');
geoshow(land,'FaceColor',[0.85 0.85 0.85],'EdgeColor','none');
plot(coastlon, coastlat, 'k-','LineWidth',0.6);

ylabel('Latitude');
xlabel('Longitude');
set(ax,'tickdir','out','FontSize',9);
xlim([20 150]); ylim([-65 0]);
apply_geo_ticks(ax, lonTicks, latTicks, true, true);    % show X and Y

text(0.83, 0.6, '$\beta S''|_{p}$','Units','normalized', ...
     'FontSize',11,'FontWeight','bold','Color','k','Interpreter','latex');

colormap(ax, cmap_modified); caxis(caxis_lim);
hold off;

% ========== 7. Panel (e): beta*S spice 300-1000m ==========
ax = subplot('Position', positions(5,:));
s_spice_bt_300_1000 = squeeze(nanmean(beta_mean(:,:,25:44),3)) .* ...
                      squeeze(nanmean(Argo.s_spice_z(:,:,25:44,:),3));
s_spice_bt_ds = deseason(fillmissing(s_spice_bt_300_1000,'linear',3), ...
                         Argo.t2, 'monthly');
%[s_spice_bt_tr, p] = trend(s_spice_bt_ds, 12, 'omitnan');
[s_spice_bt_tr, p] = calc_trend_ess(s_spice_bt_ds);

pcolor(Argo.longitude, Argo.latitude, s_spice_bt_tr'); shading flat;
hold on;
stipple(Argo.longitude, Argo.latitude, (p>=0.05)', ...
        'color','k','density',50,'marker','.','markersize',5);
[c,h] = contour(Argo.longitude, Argo.latitude, ...
                squeeze(nanmean(Argo.gm(:,:,25:44),3))', ...
                'k-','LineWidth',0.5);
clabel(c,h,'FontSize',7,'Color','k');
geoshow(land,'FaceColor',[0.85 0.85 0.85],'EdgeColor','none');
plot(coastlon, coastlat, 'k-','LineWidth',0.6);

xlabel('Longitude');
set(ax,'tickdir','out','FontSize',9);
xlim([20 150]); ylim([-65 0]);
apply_geo_ticks(ax, lonTicks, latTicks, true, false);   % show X, hide Y

text(0.83, 0.6, '$\beta S''|_{\mathrm{spice}}$','Units','normalized', ...
     'FontSize',11,'FontWeight','bold','Color','k','Interpreter','latex');

colormap(ax, cmap_modified); caxis(caxis_lim);
hold off;

% ========== 8. Panel (f): beta*S heave 300-1000m ==========
ax = subplot('Position', positions(6,:));
s_heave_bt_300_1000 = squeeze(nanmean(beta_mean(:,:,25:44),3)) .* ...
                      squeeze(nanmean(Argo.s_heave_z(:,:,25:44,:),3));
s_heave_bt_ds = deseason(fillmissing(s_heave_bt_300_1000,'linear',3), ...
                         Argo.t2, 'monthly');
%[s_heave_bt_tr, p] = trend(-s_heave_bt_ds, 12, 'omitnan');
[s_heave_bt_tr, p] = calc_trend_ess(-s_heave_bt_ds);

pcolor(Argo.longitude, Argo.latitude, s_heave_bt_tr'); shading flat;
hold on;
stipple(Argo.longitude, Argo.latitude, (p>=0.05)', ...
        'color','k','density',50,'marker','.','markersize',5);
[c,h] = contour(Argo.longitude, Argo.latitude, ...
                squeeze(nanmean(Argo.gm(:,:,25:44),3))', ...
                'k-','LineWidth',0.5);
clabel(c,h,'FontSize',7,'Color','k');
geoshow(land,'FaceColor',[0.85 0.85 0.85],'EdgeColor','none');
plot(coastlon, coastlat, 'k-','LineWidth',0.6);

c = colorbar; c.Label.String = '[kg/m^3/Yr]'; c.TickDirection = 'out'; c.Ticks = -5e-6:1e-6:5e-6;
xlabel('Longitude');
set(ax,'tickdir','out','FontSize',9);
xlim([20 150]); ylim([-65 0]);
apply_geo_ticks(ax, lonTicks, latTicks, true, false);   % show X, hide Y

text(0.83, 0.6, '$\beta S''|_{\mathrm{heave}}$','Units','normalized', ...
     'FontSize',11,'FontWeight','bold','Color','k','Interpreter','latex');

colormap(ax, cmap_modified); caxis(caxis_lim);
hold off;

% ========== 9. Add (a), (b), ... labels above each subplot ==========
for idx = 1:nrows*ncols
    x_left = positions(idx, 1);
    y_top  = positions(idx, 2) + positions(idx, 4);
    label_str = ['(' char(96+idx) ')'];
    annotation('textbox', [x_left + 0.01, y_top + 0.04, 0, 0], ...
               'String', label_str, ...
               'FitBoxToText', 'on', ...
               'EdgeColor', 'none', ...
               'FontSize', 12, ...
               'FontWeight', 'bold');
end

% ========== 10. Optional depth label ==========
% for idx = 1:6
%     subplot('Position', positions(idx,:));
%     text(0.98, 0.02, '300-1000 m','Units','normalized', ...
%          'FontSize',8,'FontStyle','italic', ...
%          'HorizontalAlignment','right','VerticalAlignment','bottom');
% end
%% Figure S1 trend of decomposed term on 1000-2000 m
%% ========== 1. Figure setup ==========
figure('Units', 'inches', 'Position', [1 1 8.5 7], 'PaperPositionMode', 'auto');

% Define margins and gaps (normalized units)
left_margin   = 0.06;
right_margin  = 0.08;
bottom_margin = 0.1;
top_margin    = 0.06;      % Space for labels above subplots
h_gap         = 0.02;     % Horizontal gap between columns
v_gap         = 0.04;     % Vertical gap between rows

nrows = 2;
ncols = 3;

sub_w = (1 - left_margin - right_margin - (ncols-1)*h_gap) / ncols;
sub_h = (1 - bottom_margin - top_margin - (nrows-1)*v_gap) / nrows;

% Pre-compute subplot positions: [x, y, width, height]
positions = zeros(nrows*ncols, 4);
for row = 1:nrows
    for col = 1:ncols
        idx = (row-1)*ncols + col;
        x = left_margin + (col-1)*(sub_w + h_gap);
        y = bottom_margin + (nrows-row)*(sub_h + v_gap);
        positions(idx, :) = [x, y, sub_w, sub_h];
    end
end

% Load land and coastline data (assumes Mapping Toolbox)
land = shaperead('landareas.shp', 'UseGeoCoords', true);
load coastlines

% ========== 2. Define common plotting parameters ==========
caxis_lim = [-1e-6 1e-6];   % Adjusted for 1000-2000 m (very weak trends)
cmap_full = cmocean('balance', 36);
cmap_modified = cmap_full([1:16, 21:36], :);  % 15 blue + 15 red, balanced

% Define standard geographic ticks for apply_geo_ticks
lon_ticks = 20:20:140;
lat_ticks = -60:20:0;

% ========== 3. Panel (a): alpha*T|z 1000-2000m ==========
ax = subplot('Position', positions(1,:));
t_z_ap_1000_2000 = squeeze(nanmean(alpha_mean(:,:,44:end),3)) .* squeeze(nanmean(Argo.t_z(:,:,44:end,:),3));
t_z_ap_ds = deseason(t_z_ap_1000_2000, Argo.t2, 'monthly');
[t_z_ap_tr, p] = calc_trend_ess(t_z_ap_ds);
pcolor(Argo.longitude, Argo.latitude, t_z_ap_tr'); shading flat;
hold on;
stipple(Argo.longitude, Argo.latitude, (p>=0.05)', 'color', 'k', 'density', 50, 'marker', '.', 'markersize', 5);
[c,h] = contour(Argo.longitude, Argo.latitude, squeeze(nanmean(Argo.gm(:,:,44:end),3))', 'k-', 'LineWidth', 0.5);
clabel(c,h, 'FontSize', 7, 'Color', 'k');
geoshow(land, 'FaceColor', [0.85 0.85 0.85], 'EdgeColor', 'none');
plot(coastlon, coastlat, 'k-', 'LineWidth', 0.6);
ylabel('Latitude');
set(gca, 'tickdir', 'out', 'FontSize', 9);
xlim([20 150]); ylim([-65 0]);
text(0.83, 0.6, '$\alpha T''|_{p}$', 'Units', 'normalized', 'FontSize', 11, 'FontWeight', 'bold', 'Color', 'k', 'Interpreter', 'latex');
colormap(gca, cmap_modified); caxis(caxis_lim);
% Apply geo ticks (hide X, show Y)
apply_geo_ticks(ax, lon_ticks, lat_ticks, false, true);
hold off;

% ========== 4. Panel (b): alpha*T spice 1000-2000m ==========
ax = subplot('Position', positions(2,:));
t_spice_ap_1000_2000 = squeeze(nanmean(alpha_mean(:,:,44:end),3)) .* squeeze(nanmean(Argo.t_spice_z(:,:,44:end,:),3));
t_spice_ap_ds = deseason(fillmissing(t_spice_ap_1000_2000, 'linear', 3), Argo.t2, 'monthly');
[t_spice_ap_tr, p] = calc_trend_ess(t_spice_ap_ds);
pcolor(Argo.longitude, Argo.latitude, t_spice_ap_tr'); shading flat;
hold on;
stipple(Argo.longitude, Argo.latitude, (p>=0.05)', 'color', 'k', 'density', 50, 'marker', '.', 'markersize', 5);
[c,h] = contour(Argo.longitude, Argo.latitude, squeeze(nanmean(Argo.gm(:,:,44:end),3))', 'k-', 'LineWidth', 0.5);
clabel(c,h, 'FontSize', 7, 'Color', 'k');
geoshow(land, 'FaceColor', [0.85 0.85 0.85], 'EdgeColor', 'none');
plot(coastlon, coastlat, 'k-', 'LineWidth', 0.6);
set(gca, 'tickdir', 'out', 'FontSize', 9);
xlim([20 150]); ylim([-65 0]);
text(0.83, 0.6, '$\alpha T''|_{\mathrm{spice}}$', 'Units', 'normalized', 'FontSize', 11, 'FontWeight', 'bold', 'Color', 'k', 'Interpreter', 'latex');
colormap(gca, cmap_modified); caxis(caxis_lim);
% Apply geo ticks (hide X, hide Y)
apply_geo_ticks(ax, lon_ticks, lat_ticks, false, false);
hold off;

% ========== 5. Panel (c): alpha*T heave 1000-2000m ==========
ax = subplot('Position', positions(3,:));
t_heave_ap_1000_2000 = squeeze(nanmean(alpha_mean(:,:,44:end),3)) .* squeeze(nanmean(Argo.t_heave_z(:,:,44:end,:),3));
t_heave_ap_ds = deseason(fillmissing(t_heave_ap_1000_2000, 'linear', 3), Argo.t2, 'monthly');
[t_heave_ap_tr, p] = calc_trend_ess(t_heave_ap_ds);
pcolor(Argo.longitude, Argo.latitude, -t_heave_ap_tr'); shading flat;
hold on;
stipple(Argo.longitude, Argo.latitude, (p>=0.05)', 'color', 'k', 'density', 50, 'marker', '.', 'markersize', 5);
[c,h] = contour(Argo.longitude, Argo.latitude, squeeze(nanmean(Argo.gm(:,:,44:end),3))', 'k-', 'LineWidth', 0.5);
clabel(c,h, 'FontSize', 7, 'Color', 'k');
geoshow(land, 'FaceColor', [0.85 0.85 0.85], 'EdgeColor', 'none');
plot(coastlon, coastlat, 'k-', 'LineWidth', 0.6);
c = colorbar;
c.Label.String = '[kg/m^3/Yr]';
c.TickDirection = 'out'; 
set(gca, 'tickdir', 'out', 'FontSize', 9);
xlim([20 150]); ylim([-65 0]);
text(0.83, 0.6, '$\alpha T''|_{\mathrm{heave}}$', 'Units', 'normalized', 'FontSize', 11, 'FontWeight', 'bold', 'Color', 'k', 'Interpreter', 'latex');
colormap(gca, cmap_modified); caxis(caxis_lim);
% Apply geo ticks (hide X, hide Y)
apply_geo_ticks(ax, lon_ticks, lat_ticks, false, false);
hold off;

% ========== 6. Panel (d): beta*S|z 1000-2000m ==========
ax = subplot('Position', positions(4,:));
s_z_bt_1000_2000 = squeeze(nanmean(beta_mean(:,:,44:end),3)) .* squeeze(nanmean(Argo.s_z(:,:,44:end,:),3));
s_z_bt_ds = deseason(s_z_bt_1000_2000, Argo.t2, 'monthly');
[s_z_bt_tr, p] = calc_trend_ess(s_z_bt_ds);
pcolor(Argo.longitude, Argo.latitude, s_z_bt_tr'); shading flat;
hold on;
stipple(Argo.longitude, Argo.latitude, (p>=0.05)', 'color', 'k', 'density', 50, 'marker', '.', 'markersize', 5);
[c,h] = contour(Argo.longitude, Argo.latitude, squeeze(nanmean(Argo.gm(:,:,44:end),3))', 'k-', 'LineWidth', 0.5);
clabel(c,h, 'FontSize', 7, 'Color', 'k');
geoshow(land, 'FaceColor', [0.85 0.85 0.85], 'EdgeColor', 'none');
plot(coastlon, coastlat, 'k-', 'LineWidth', 0.6);
ylabel('Latitude');
xlabel('Longitude');
set(gca, 'tickdir', 'out', 'FontSize', 9);
xlim([20 150]); ylim([-65 0]);
text(0.83, 0.6, '$\beta S''|_{p}$', 'Units', 'normalized', 'FontSize', 11, 'FontWeight', 'bold', 'Color', 'k', 'Interpreter', 'latex');
colormap(gca, cmap_modified); caxis(caxis_lim);
% Apply geo ticks (show X, show Y)
apply_geo_ticks(ax, lon_ticks, lat_ticks, true, true);
hold off;

% ========== 7. Panel (e): beta*S spice 1000-2000m ==========
ax = subplot('Position', positions(5,:));
s_spice_bt_1000_2000 = squeeze(nanmean(beta_mean(:,:,44:end),3)) .* squeeze(nanmean(Argo.s_spice_z(:,:,44:end,:),3));
s_spice_bt_ds = deseason(fillmissing(s_spice_bt_1000_2000, 'linear', 3), Argo.t2, 'monthly');
[s_spice_bt_tr, p] = calc_trend_ess(s_spice_bt_ds);
pcolor(Argo.longitude, Argo.latitude, s_spice_bt_tr'); shading flat;
hold on;
stipple(Argo.longitude, Argo.latitude, (p>=0.05)', 'color', 'k', 'density', 50, 'marker', '.', 'markersize', 5);
[c,h] = contour(Argo.longitude, Argo.latitude, squeeze(nanmean(Argo.gm(:,:,44:end),3))', 'k-', 'LineWidth', 0.5);
clabel(c,h, 'FontSize', 7, 'Color', 'k');
geoshow(land, 'FaceColor', [0.85 0.85 0.85], 'EdgeColor', 'none');
plot(coastlon, coastlat, 'k-', 'LineWidth', 0.6);
xlabel('Longitude');
set(gca, 'tickdir', 'out', 'FontSize', 9);
xlim([20 150]); ylim([-65 0]);
text(0.83, 0.6, '$\beta S''|_{\mathrm{spice}}$', 'Units', 'normalized', 'FontSize', 11, 'FontWeight', 'bold', 'Color', 'k', 'Interpreter', 'latex');
colormap(gca, cmap_modified); caxis(caxis_lim);
% Apply geo ticks (show X, hide Y)
apply_geo_ticks(ax, lon_ticks, lat_ticks, true, false);
hold off;

% ========== 8. Panel (f): beta*S heave 1000-2000m ==========
ax = subplot('Position', positions(6,:));
s_heave_bt_1000_2000 = squeeze(nanmean(beta_mean(:,:,44:end),3)) .* squeeze(nanmean(Argo.s_heave_z(:,:,44:end,:),3));
s_heave_bt_ds = deseason(fillmissing(s_heave_bt_1000_2000, 'linear', 3), Argo.t2, 'monthly');
[s_heave_bt_tr, p] = calc_trend_ess(-s_heave_bt_ds);
pcolor(Argo.longitude, Argo.latitude, s_heave_bt_tr'); shading flat;
hold on;
stipple(Argo.longitude, Argo.latitude, (p>=0.05)', 'color', 'k', 'density', 50, 'marker', '.', 'markersize', 5);
[c,h] = contour(Argo.longitude, Argo.latitude, squeeze(nanmean(Argo.gm(:,:,44:end),3))', 'k-', 'LineWidth', 0.5);
clabel(c,h, 'FontSize', 7, 'Color', 'k');
geoshow(land, 'FaceColor', [0.85 0.85 0.85], 'EdgeColor', 'none');
plot(coastlon, coastlat, 'k-', 'LineWidth', 0.6);
c = colorbar;
c.Label.String = '[kg/m^3/Yr]';
c.TickDirection = 'out'; 
xlabel('Longitude');
set(gca, 'tickdir', 'out', 'FontSize', 9);
xlim([20 150]); ylim([-65 0]);
text(0.83, 0.6, '$\beta S''|_{\mathrm{heave}}$', 'Units', 'normalized', 'FontSize', 11, 'FontWeight', 'bold', 'Color', 'k', 'Interpreter', 'latex');
colormap(gca, cmap_modified); caxis(caxis_lim);
% Apply geo ticks (show X, hide Y)
apply_geo_ticks(ax, lon_ticks, lat_ticks, true, false);
hold off;

% ========== 9. Add labels (a), (b), (c), ... ABOVE each subplot ==========
for idx = 1:nrows*ncols
    x_left = positions(idx, 1);
    y_top  = positions(idx, 2) + positions(idx, 4);
    label_str = ['(' char(96+idx) ')'];  % a=97, b=98, ...
    annotation('textbox', [x_left + 0.01, y_top + 0.04, 0, 0], ...
               'String', label_str, ...
               'FitBoxToText', 'on', ...
               'EdgeColor', 'none', ...
               'FontSize', 12, ...
               'FontWeight', 'bold');
end

%% plot figure 4: verticle trend:
%% prepare data:
% Initialize matrices to store alpha and beta values
[nlon, nlat, ndepth, ntime] = size(Argo.CT_month_all);
alpha = nan(nlon, nlat, ndepth, ntime);
beta = nan(nlon, nlat, ndepth, ntime);

% Loop through time and depth
for t = 1:ntime
    for d = 1:ndepth
        % Extract sea pressure at this depth level (constant across lat, lon)
        p_level = Argo.pressure(d); % sea pressure for the specific depth in decibars
        
        % Compute alpha (thermal expansion coefficient)
        alpha(:, :, d, t) = gsw_alpha(Argo.SA_month_all(:, :, d, t), Argo.CT_month_all(:, :, d, t), p_level);
        
        % Compute beta (saline contraction coefficient)
        beta(:, :, d, t) = gsw_beta(Argo.SA_month_all(:, :, d, t), Argo.CT_month_all(:, :, d, t), p_level);
    end
end

%prepare trends:
nlon   = size(Argo.t_z,1);
nlat   = size(Argo.t_z,2);
ndepth = size(Argo.t_z,3);
ntime  = size(Argo.t_z,4);

% Preallocate trend arrays: depth is 3rd dimension
t_z_tr        = nan(nlon,nlat,ndepth);
t_z_p         = nan(nlon,nlat,ndepth);
t_spice_z_tr  = nan(nlon,nlat,ndepth);
t_spice_z_p   = nan(nlon,nlat,ndepth);
t_heave_z_tr  = nan(nlon,nlat,ndepth);
t_heave_z_p   = nan(nlon,nlat,ndepth);
s_z_tr        = nan(nlon,nlat,ndepth);
s_z_p         = nan(nlon,nlat,ndepth);
s_spice_z_tr  = nan(nlon,nlat,ndepth);
s_spice_z_p   = nan(nlon,nlat,ndepth);
s_heave_z_tr  = nan(nlon,nlat,ndepth);
s_heave_z_p   = nan(nlon,nlat,ndepth);

%%
for i_lo = 1:nlon
    for i_la = 1:nlat
        for i_dep = 1:ndepth
            % Extract time series at this depth
            temp_series      = squeeze(Argo.t_z(i_lo,i_la,i_dep,:));
            temp_spice       = squeeze(Argo.t_spice_z(i_lo,i_la,i_dep,:));
            temp_heave       = squeeze(Argo.t_heave_z(i_lo,i_la,i_dep,:));
            salt_series      = squeeze(Argo.s_z(i_lo,i_la,i_dep,:));
            salt_spice       = squeeze(Argo.s_spice_z(i_lo,i_la,i_dep,:));
            salt_heave       = squeeze(Argo.s_heave_z(i_lo,i_la,i_dep,:));

            % Fill missing values
            temp_series      = fillmissing(temp_series,'linear');
            temp_spice       = fillmissing(temp_spice,'linear');
            temp_heave       = fillmissing(temp_heave,'linear');
            salt_series      = fillmissing(salt_series,'linear');
            salt_spice       = fillmissing(salt_spice,'linear');
            salt_heave       = fillmissing(salt_heave,'linear');

            % Deseason
            temp_series      = deseason(temp_series,Argo.t2,12);
            temp_spice       = deseason(temp_spice,Argo.t2,12);
            temp_heave       = deseason(temp_heave,Argo.t2,12);
            salt_series      = deseason(salt_series,Argo.t2,12);
            salt_spice       = deseason(salt_spice,Argo.t2,12);
            salt_heave       = deseason(salt_heave,Argo.t2,12);

            % Compute trend along time (output scalar per depth)
            % [t_z_tr(i_lo,i_la,i_dep), t_z_p(i_lo,i_la,i_dep)]               = trend(temp_series,12,'omitnan');
            % [t_spice_z_tr(i_lo,i_la,i_dep), t_spice_z_p(i_lo,i_la,i_dep)]   = trend(temp_spice,12,'omitnan');
            % [t_heave_z_tr(i_lo,i_la,i_dep), t_heave_z_p(i_lo,i_la,i_dep)]   = trend(temp_heave,12,'omitnan');
            % [s_z_tr(i_lo,i_la,i_dep), s_z_p(i_lo,i_la,i_dep)]               = trend(salt_series,12,'omitnan');
            % [s_spice_z_tr(i_lo,i_la,i_dep), s_spice_z_p(i_lo,i_la,i_dep)]   = trend(salt_spice,12,'omitnan');
            % [s_heave_z_tr(i_lo,i_la,i_dep), s_heave_z_p(i_lo,i_la,i_dep)]   = trend(salt_heave,12,'omitnan');

            [t_z_tr(i_lo,i_la,i_dep), t_z_p(i_lo,i_la,i_dep)]               = calc_trend_ess_1d(temp_series);
            [t_spice_z_tr(i_lo,i_la,i_dep), t_spice_z_p(i_lo,i_la,i_dep)]   = calc_trend_ess_1d(temp_spice);
            [t_heave_z_tr(i_lo,i_la,i_dep), t_heave_z_p(i_lo,i_la,i_dep)]   = calc_trend_ess_1d(temp_heave);
            [s_z_tr(i_lo,i_la,i_dep), s_z_p(i_lo,i_la,i_dep)]               = calc_trend_ess_1d(salt_series);
            [s_spice_z_tr(i_lo,i_la,i_dep), s_spice_z_p(i_lo,i_la,i_dep)]   = calc_trend_ess_1d(salt_spice);
            [s_heave_z_tr(i_lo,i_la,i_dep), s_heave_z_p(i_lo,i_la,i_dep)]   = calc_trend_ess_1d(salt_heave);
        end
    end
end

beta_mean=nanmean(beta,4);
alpha_mean=nanmean(alpha,4);
Argo.sal=nanmean(Argo.SA_month_all,4);
%% ========== Figure setup ========== (Figure 4)
figure('Position', [100 100 650 750]);

% ========== Define subplot positions ==========
x_positions = [0.08, 0.37, 0.67];
width = 0.27;
height = 0.18;

esio_y_top = 0.75;
esio_y_bot = 0.52;
wsio_y_top = 0.29;
wsio_y_bot = 0.08;

block_y = [esio_y_top, esio_y_bot;   % ESIO
           wsio_y_top, wsio_y_bot];  % WSIO

subplot_positions = zeros(12, 4);
for block = 1:2
    for row = 1:2
        for col = 1:3
            idx = (block-1)*6 + (row-1)*3 + col;
            subplot_positions(idx, :) = [x_positions(col), block_y(block, row), width, height];
        end
    end
end

panel_labels = {'(a)', '(b)', '(c)', '(d)', '(e)', '(f)', ...
                '(g)', '(h)', '(i)', '(j)', '(k)', '(l)'};

% ========== Plotting parameters ==========
caxis_lim = [-6e-6, 6e-6];
density_contours = [25, 26.2, 26.8, 27:0.1:28];
salinity_contours = [34.8, 35.8];
stipple_density = 50;

% Latitude ticks
lat_fine = -65:5:-5;

% ========== CREATE GRIDS FOR STIPPLE (FIX) ==========
[lat_grid, pres_grid] = meshgrid(Argo.latitude, Argo.pressure);

% ========== Region loop ==========
regions = {[90, 120], [40, 70]};  % ESIO, WSIO

for r = 1:2
    lon_lim = regions{r};
    i_lo = find(Argo.longitude >= lon_lim(1) & Argo.longitude <= lon_lim(2));
    
    alpha_mean_region = squeeze(nanmean(alpha_mean(i_lo, :, :), 1));
    beta_mean_region  = squeeze(nanmean(beta_mean(i_lo, :, :), 1));
    
    for row = 1:2
        for col = 1:3
            sub_idx = (r-1)*6 + (row-1)*3 + col;
            subplot('Position', subplot_positions(sub_idx, :));
            
            % --- Select data and set titles ---
            if row == 1  % Temperature
                if col == 1
                    data = squeeze(nanmean(t_z_tr(i_lo, :, :), 1))' .* alpha_mean_region';
                    sig  = squeeze(nanmean(t_z_p(i_lo, :, :), 1))' >= 0.05;
                    title_str = '\alpha T''|p'; 
                elseif col == 2
                    data = squeeze(nanmean(t_spice_z_tr(i_lo, :, :), 1))' .* alpha_mean_region';
                    sig  = squeeze(nanmean(t_spice_z_p(i_lo, :, :), 1))' >= 0.05;
                    title_str = '\alpha T''|Spice';
                else
                    data = -squeeze(nanmean(t_heave_z_tr(i_lo, :, :), 1))' .* alpha_mean_region';
                    sig  = squeeze(nanmean(t_heave_z_p(i_lo, :, :), 1))' >= 0.05;
                    title_str = '\alpha T''|Heave'; 
                end
            else  % Salinity
                if col == 1
                    data = squeeze(nanmean(s_z_tr(i_lo, :, :), 1))' .* beta_mean_region';
                    sig  = squeeze(nanmean(s_z_p(i_lo, :, :), 1))' >= 0.05;
                    title_str = '\beta S''|p'; 
                elseif col == 2
                    data = squeeze(nanmean(s_spice_z_tr(i_lo, :, :), 1))' .* beta_mean_region';
                    sig  = squeeze(nanmean(s_spice_z_p(i_lo, :, :), 1))' >= 0.05;
                    title_str = '\beta S''|Spice'; 
                else
                    data = -squeeze(nanmean(s_heave_z_tr(i_lo, :, :), 1))' .* beta_mean_region';
                    sig  = squeeze(nanmean(s_heave_z_p(i_lo, :, :), 1))' >= 0.05;
                    title_str = '\beta S''|Heave'; 
                end
            end
            
            % --- Plot ---
            pcolor(Argo.latitude, Argo.pressure, data);
            shading flat;
            hold on;
            
            % --- Stipple (use grids) ---
            stipple(lat_grid, pres_grid, sig, ...
                    'color', [0.4,0.4,0.4], 'density', stipple_density, ...
                    'marker', '.', 'markersize', 4);
            
            % --- Density contours ---
            [c, h] = contour(Argo.latitude, Argo.pressure, ...
                squeeze(nanmean(Argo.gm(i_lo, :, :), 1))', ...
                density_contours, 'w');
            clabel(c, h, 'Color', 'w', 'LabelSpacing', 200);
            
            % --- Salinity contours ---
            [cs, hs] = contour(Argo.latitude, Argo.pressure, ...
                squeeze(nanmean(Argo.sal(i_lo, :, :), 1))', ...
                salinity_contours, 'k', 'LineWidth', 1);
            clabel(cs, hs, 'Color', 'k', 'LabelSpacing', 200);
            
            % --- Axis labels ---
            if col == 1
                ylabel('Pressure (dbar)');
            else
                ylabel('');
            end
            
            if (r == 2 && row == 2)
                xlabel('Latitude');
            else
                xlabel('');
            end
            
            % --- Panel label ---
            text(Argo.latitude(1), -128, panel_labels{sub_idx}, ...
                 'FontWeight', 'bold', 'FontSize', 12);
            
            % --- Formatting ---
            set(gca, 'ydir', 'reverse', 'tickdir', 'out', 'FontSize', 9);
            xlim([-65 -5]);
            ylim([0 2000]);
            
            % --- Tick marks: custom °S labels only on bottom row ---
            set(gca, 'XTick', lat_fine);
            if (r == 2 && row == 2)
                % Convert negative numbers to absolute values with °S
                lat_labels = arrayfun(@(x) sprintf('%d°S', abs(x)), lat_fine, 'UniformOutput', false);
                set(gca, 'XTickLabel', lat_labels);
            else
                set(gca, 'XTickLabel', {});
            end
            
            % Hide Y-axis tick labels for non-first columns
            if col ~= 1
                set(gca,'YTickLabel',{});
            end
            
            title(title_str);
            cmocean('balance', 30);
            caxis(caxis_lim);
            
            hold off;
        end
    end
end

% ========== Shared colourbar ==========
h = colorbar('Position', [0.95 0.06 0.015 0.88]);
h.TickDirection = 'out'; % Colorbar ticks point outward
ylabel(h, 'kg m^{-3} yr^{-1}');

cmap = cmocean('balance', 32);
cmap = cmap([1:15, 18:32], :);
colormap(cmap);


%% Figure for paper2 figure 5 EOF ESIO
%prepare the data
%% calculate mean terms for adding contours
%ESIO
 i_lo=find(Argo.longitude >=90 & Argo.longitude <=120);
SA_mean_E=squeeze(nanmean(nanmean(Argo.SA_month_all(i_lo,:,:,:),1),4));
PV_mean_E=squeeze(nanmean(nanmean(Argo.interp_PV(i_lo,:,:,:),1),4));
GM_mean_E=squeeze(nanmean(nanmean(Argo.gamma_n(i_lo,:,:,:),1),4));

i_lo=find(Argo.longitude >=40 & Argo.longitude <=70);
SA_mean_w=squeeze(nanmean(nanmean(Argo.SA_month_all(i_lo,:,:,:),1),4));
PV_mean_w=squeeze(nanmean(nanmean(Argo.interp_PV(i_lo,:,:,:),1),4));
GM_mean_w=squeeze(nanmean(nanmean(Argo.gamma_n(i_lo,:,:,:),1),4));

% Initialize matrices to store alpha and beta values
[nlon, nlat, ndepth, ntime] = size(Argo.CT_month_all);
alpha = nan(nlon, nlat, ndepth, ntime);
beta = nan(nlon, nlat, ndepth, ntime);

% Loop through time and depth
for t = 1:ntime
    for d = 1:ndepth
        % Extract sea pressure at this depth level (constant across lat, lon)
        p_level = Argo.pressure(d); % sea pressure for the specific depth in decibars
        
        % Compute alpha (thermal expansion coefficient)
        alpha(:, :, d, t) = gsw_alpha(Argo.SA_month_all(:, :, d, t), Argo.CT_month_all(:, :, d, t), p_level);
        
        % Compute beta (saline contraction coefficient)
        beta(:, :, d, t) = gsw_beta(Argo.SA_month_all(:, :, d, t), Argo.CT_month_all(:, :, d, t), p_level);
    end
end


time=Argo.t2;
load DMI.mat
load nino34.mat
%% ========== Colour definitions (your own) ==========
pc1_color = 'k';            % PC1: black solid
dmi_color = [0.85 0.6 0.1]; % DMI: yellow-orange (solid)
enso_color = 'b';           % ENSO: blue dashed

% ESIO no smoothed EOF
lo_range=[90 120];
lat=Argo.latitude;
la_range=[-65 0];
i_lat=find(Argo.latitude >= -65 & Argo.latitude <=0);
i_lo=find(Argo.longitude >=lo_range(1) & Argo.longitude <=lo_range(2));
i_d=find(Argo.pressure >= 0);
beta_sa=beta.*Argo.SA_month_all;
sa_mean=squeeze(nanmean(beta_sa(i_lo,i_lat,i_d,:),1));
%sa_mean=movmean(sa_mean,12,3);
%sa_mean=detrend(sa_mean);
sa_mean=movmean(deseason(sa_mean,Argo.t2),12,3);
%sa_mean=deseason(sa_mean,Argo.t2);
%sa_mean=detrend(sa_mean);

alpha_ct=alpha.*Argo.CT_month_all;
ct_mean=squeeze(nanmean(alpha_ct(i_lo,i_lat,i_d,:),1));
%ct_mean=movmean(ct_mean,12,3);
%ct_mean=detrend(ct_mean);
 ct_mean=movmean(deseason(ct_mean,Argo.t2),12,3);
%ct_mean=detrend(ct_mean);
%ct_mean=deseason(ct_mean,Argo.t2);

[ny, nz, nt] = size(sa_mean);
sa_mean_reshape = reshape(sa_mean, ny*nz, nt);
ct_mean_reshape = reshape(ct_mean, ny*nz, nt); 
ct_sa_reshape=[sa_mean_reshape; ct_mean_reshape];
valid = all(~isnan(ct_sa_reshape), 2);
ct_sa_reshape_valid = ct_sa_reshape(valid, :);

ct_sa_anom = detrend(ct_sa_reshape_valid', 'linear');
[EOFs, PCs, latent] = pca(ct_sa_anom);

EOF_maps = nan(ny*nz, size(EOFs, 2));
EOF_maps(valid, :) = EOFs;
EOF_maps_reshaped = reshape(EOF_maps, ny, nz*2, []);

% plot total eof on depth level T|z S|z
figure('Position',[300,50,800,700])

% ========== Row 1: Total on Pressure Levels (a-c) ==========
% (a) Temperature spatial pattern
subplot('Position',[0.07,0.7,0.2,0.26]);
pcolor(lat(i_lat), Argo.pressure(i_d), EOF_maps_reshaped(:,59:end,1)'.*std(PCs(:,1)));
shading interp;
cmap_full = cmocean('balance', 36);
cmap_full = cmap_full([1:17, 20:36], :);   
colormap(gca, cmap_full);
caxis([-5e-5 5e-5])
set(gca, 'ydir', 'reverse', 'XAxisLocation', 'bottom', 'tickdir', 'out')
hold on
[c,h]=contour(Argo.latitude,Argo.pressure,SA_mean_E',[34.2,34.7,34.8,35.8],'k');
clabel(c,h);
[c,h]=contour(Argo.latitude,Argo.pressure,GM_mean_E',[25,25.8,26,26.2,26.8,27:0.1:28],'w');
clabel(c,h,'color','w');
ylabel('Pressure (dbar)');
xlabel('');
set(gca, 'XTickLabel', {});
% Label at top-left using text (no warning, perfect placement)
text(0.02, 1.13, '(a) $T''|_p$', 'Units', 'normalized', ...
     'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', ...
     'FontSize', 11, 'FontWeight', 'bold', 'Interpreter', 'latex');

% (b) Salinity spatial pattern
subplot('Position',[0.3,0.7,0.2,0.26]);
pcolor(lat(i_lat), Argo.pressure(i_d), -EOF_maps_reshaped(:,1:58,1)'.*std(PCs(:,1)));
shading interp;
colorbar;
cmap_full = cmocean('balance', 36);
cmap_full = cmap_full([1:17, 20:36], :);   
colormap(gca, cmap_full);
caxis([-5e-5 5e-5])
set(gca, 'ydir', 'reverse', 'XAxisLocation', 'bottom', 'tickdir', 'out')
hold on
[c,h]=contour(Argo.latitude,Argo.pressure,SA_mean_E',[34.2,34.7,34.8,35.8],'k');
clabel(c,h);
[c,h]=contour(Argo.latitude,Argo.pressure,GM_mean_E',[25,25.8,26,26.2,26.8,27:0.1:28],'w');
clabel(c,h,'color','w');
ylabel('');
xlabel('');
set(gca, 'XTickLabel', {});
set(gca, 'YTickLabel', {});
text(0.02, 1.13, '(b) $S''|_p$', 'Units', 'normalized', ...
     'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', ...
     'FontSize', 11, 'FontWeight', 'bold', 'Interpreter', 'latex');

% (c) PC time series
subplot('Position',[0.58,0.75,0.4,0.18]);
plot(time, -PCs(:,1)./std(PCs(:,1)), [pc1_color '-'], 'LineWidth', 1.5);
hold on
text(0.02, 1.13, '(c) $PC1|_p$', 'Units', 'normalized', ...
     'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', ...
     'FontSize', 11, 'FontWeight', 'bold', 'Interpreter', 'latex');
subtitle(['Explained variance: ',num2str(round(latent(1)/sum(latent)*100,1)),'%'])
xlabel('Time');
ylabel('Amplitude');
ylim([-2.5 2.5])
set(gca, 'tickdir', 'out')

pc1=PCs(:,1)./std(PCs(:,1));
pc1=(pc1-nanmean(pc1))./std(pc1);

% DMI: your own colour (yellow-orange)
dmi_smooth = movmean((DMI.Index_Argo-nanmean(DMI.Index_Argo))./std(DMI.Index_Argo,'omitnan'),12,1);
plot(DMI.Dates_Argo, dmi_smooth, 'Color', dmi_color, 'LineWidth', 1.5);
correlation_dmi = corr(pc1, dmi_smooth, 'rows', 'complete');
text(DMI.Dates_Argo(10), -2.2, ['Corr with DMI: ', num2str(round(correlation_dmi,2))]);

% ENSO: blue dashed
enso_smooth = movmean((nino34.Index_Argo-nanmean(nino34.Index_Argo))./std(nino34.Index_Argo,'omitnan'),12,1);
plot(nino34.Dates_Argo, enso_smooth, [enso_color '--'], 'LineWidth', 1.5);
correlation_enso = corr(pc1, enso_smooth, 'rows', 'complete');
text(nino34.Dates_Argo(120), -2.2, ['Corr with ENSO: ', num2str(round(correlation_enso,2))]);

% Grey zero line
plot(time, zeros(size(time)), 'Color', [0.5 0.5 0.5], 'LineStyle', '--', 'LineWidth', 0.8);
legend('PC1','DMI','ENSO', 'Location', 'best');
hold off


% ========== Row 2: Spice EOF (d-f) ==========
lo_range=[90 120];
lat=Argo.latitude;
la_range=[-65 0];
i_lat=find(Argo.latitude >= la_range(1) & Argo.latitude <=la_range(2));
i_lo=find(Argo.longitude >=lo_range(1) & Argo.longitude <=lo_range(2));
beta_sa=beta.*Argo.s_spice_z;
sa_mean=squeeze(nanmean(beta_sa(i_lo,i_lat,:,:),1));
%sa_mean=movmean(sa_mean,12,3);
%sa_mean=detrend(sa_mean);
sa_mean=movmean(deseason(sa_mean,Argo.t2),12,3);
%sa_mean=deseason(sa_mean,Argo.t2);
%sa_mean=detrend(sa_mean);

alpha_ct=alpha.*Argo.t_spice_z;
ct_mean=squeeze(nanmean(alpha_ct(i_lo,i_lat,:,:),1));
%ct_mean=movmean(ct_mean,12,3);
%ct_mean=detrend(ct_mean);
ct_mean=movmean(deseason(ct_mean,Argo.t2),12,3);
%ct_mean=detrend(ct_mean);
%ct_mean=deseason(ct_mean,Argo.t2);

[ny, nz, nt] = size(sa_mean);
sa_mean_reshape = reshape(sa_mean, ny*nz, nt);
valid_sa = all(~isnan(sa_mean_reshape), 2);
sa_reshape_valid=sa_mean_reshape(valid_sa,:);
ct_mean_reshape = reshape(ct_mean, ny*nz, nt);
valid_ct = all(~isnan(ct_mean_reshape), 2);
ct_reshape_valid=ct_mean_reshape(valid_ct,:);
ct_sa_reshape=[sa_reshape_valid; ct_reshape_valid];

ct_sa_anom = detrend(ct_sa_reshape', 'linear');
[EOFs, PCs, latent] = pca(ct_sa_anom);

EOF_maps = nan(ny*nz*2, size(EOFs, 2));
EOF_maps([valid_ct,valid_sa], :) = EOFs;
EOF_maps_reshaped = reshape(EOF_maps, ny, nz*2, []);

% (d) Temperature pattern (Spice)
subplot('Position',[0.07,0.38,0.2,0.26]);
pcolor(lat(i_lat), Argo.pressure(i_d), -EOF_maps_reshaped(:,59:end,1)'.*std(PCs(:,1)));
shading interp;
cmap_full = cmocean('balance', 36);
cmap_full = cmap_full([1:17, 20:36], :);   
colormap(gca, cmap_full);
caxis([-5e-5 5e-5])
set(gca, 'ydir', 'reverse', 'XAxisLocation', 'bottom', 'tickdir', 'out')
hold on
[c,h]=contour(Argo.latitude,Argo.pressure,SA_mean_E',[34.2,34.7,34.8,35.8],'k');
clabel(c,h);
[c,h]=contour(Argo.latitude,Argo.pressure,GM_mean_E',[25,25.8,26,26.2,26.8,27:0.1:28],'w');
clabel(c,h,'color','w');
ylabel('Pressure (dbar)');
xlabel('');
set(gca, 'XTickLabel', {});
text(0.02, 1.13, '(d) $T''|_{\rm spice}$', 'Units', 'normalized', ...
     'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', ...
     'FontSize', 11, 'FontWeight', 'bold', 'Interpreter', 'latex');

% (e) Salinity pattern (Spice)
subplot('Position',[0.3,0.38,0.2,0.26]);
pcolor(lat(i_lat), Argo.pressure(i_d), -EOF_maps_reshaped(:,1:58,1)'.*std(PCs(:,1)));
shading interp;
colorbar;
cmap_full = cmocean('balance', 36);
cmap_full = cmap_full([1:17, 20:36], :);   
colormap(gca, cmap_full);
caxis([-5e-5 5e-5])
set(gca, 'ydir', 'reverse', 'XAxisLocation', 'bottom', 'tickdir', 'out')
hold on
[c,h]=contour(Argo.latitude,Argo.pressure,SA_mean_E',[34.2,34.7,34.8,35.8],'k');
clabel(c,h);
[c,h]=contour(Argo.latitude,Argo.pressure,GM_mean_E',[25,25.8,26,26.2,26.8,27:0.1:28],'w');
clabel(c,h,'color','w');
ylabel('');
xlabel('');
set(gca, 'XTickLabel', {});
set(gca, 'YTickLabel', {});
text(0.02, 1.13, '(e) $S''|_{\rm spice}$', 'Units', 'normalized', ...
     'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', ...
     'FontSize', 11, 'FontWeight', 'bold', 'Interpreter', 'latex');

% (f) PC time series (Spice)
subplot('Position',[0.58,0.4,0.4,0.18]);
plot(time, -PCs(:,1)./std(PCs(:,1)), [pc1_color '-'], 'LineWidth', 1.5);
hold on
text(0.02, 1.13, '(f) $PC1|_{\rm spice}$', 'Units', 'normalized', ...
     'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', ...
     'FontSize', 11, 'FontWeight', 'bold', 'Interpreter', 'latex');
subtitle(['Explained variance: ',num2str(round(latent(1)/sum(latent)*100,1)),'%'])
xlabel('Time');
ylabel('Amplitude');
ylim([-2.5 2.5])
set(gca, 'tickdir', 'out')

pc1=PCs(:,1)./std(PCs(:,1));
pc1=(pc1-nanmean(pc1))./std(pc1);

dmi_smooth = movmean((DMI.Index_Argo-nanmean(DMI.Index_Argo))./std(DMI.Index_Argo,'omitnan'),12,1);
plot(DMI.Dates_Argo, dmi_smooth, 'Color', dmi_color, 'LineWidth', 1.5);
correlation_dmi = corr(pc1, dmi_smooth, 'rows', 'complete');
text(DMI.Dates_Argo(10), -2.2, ['Corr with DMI: ', num2str(round(correlation_dmi,2))]);

enso_smooth = movmean((nino34.Index_Argo-nanmean(nino34.Index_Argo))./std(nino34.Index_Argo,'omitnan'),12,1);
plot(nino34.Dates_Argo, enso_smooth, [enso_color '--'], 'LineWidth', 1.5);
correlation_enso = corr(pc1, enso_smooth, 'rows', 'complete');
text(nino34.Dates_Argo(120), -2.2, ['Corr with ENSO: ', num2str(round(correlation_enso,2))]);

plot(time, zeros(size(time)), 'Color', [0.5 0.5 0.5], 'LineStyle', '--', 'LineWidth', 0.8);
legend('PC1','DMI','ENSO', 'Location', 'best');
hold off


% ========== Row 3: Heave EOF (g-i) ==========
lo_range=[90 120];
lat=Argo.latitude;
la_range=[-65 0];
i_lat=find(Argo.latitude >= la_range(1) & Argo.latitude <=la_range(2));
i_lo=find(Argo.longitude >=lo_range(1) & Argo.longitude <=lo_range(2));
beta_sa=beta.*Argo.s_heave_z;
sa_mean=squeeze(nanmean(beta_sa(i_lo,i_lat,:,:),1));
%sa_mean=movmean(sa_mean,12,3);
%sa_mean=detrend(sa_mean);
 sa_mean=movmean(deseason(sa_mean,Argo.t2),12,3);
%sa_mean=detrend(sa_mean,'linear');
%sa_mean=deseason(sa_mean,Argo.t2);

alpha_ct=alpha.*Argo.t_heave_z;
ct_mean=squeeze(nanmean(alpha_ct(i_lo,i_lat,:,:),1));
%ct_mean=movmean(ct_mean,12,3);
%ct_mean=detrend(ct_mean);
ct_mean=movmean(deseason(ct_mean,Argo.t2),12,3);
%ct_mean=detrend(ct_mean,'linear');
%ct_mean=deseason(ct_mean,Argo.t2);

[ny, nz, nt] = size(sa_mean);
sa_mean_reshape = reshape(sa_mean, ny*nz, nt);
valid_sa = all(~isnan(sa_mean_reshape), 2);
sa_reshape_valid=sa_mean_reshape(valid_sa,:);
ct_mean_reshape = reshape(ct_mean, ny*nz, nt);
valid_ct = all(~isnan(ct_mean_reshape), 2);
ct_reshape_valid=ct_mean_reshape(valid_ct,:);
ct_sa_reshape=[sa_reshape_valid; ct_reshape_valid];

ct_sa_anom = detrend(ct_sa_reshape', 'linear');
[EOFs, PCs, latent] = pca(ct_sa_anom);

EOF_maps = nan(ny*nz*2, size(EOFs, 2));
EOF_maps([valid_ct,valid_sa], :) = EOFs;
EOF_maps_reshaped = reshape(EOF_maps, ny, nz*2, []);

% (g) Temperature pattern (Heave)
subplot('Position',[0.07,0.07,0.2,0.26]);
pcolor(lat(i_lat), Argo.pressure(i_d), EOF_maps_reshaped(:,59:end,1)'.*std(PCs(:,1)));
shading interp;
cmap_full = cmocean('balance', 36);
cmap_full = cmap_full([1:17, 20:36], :);   
colormap(gca, cmap_full);
caxis([-5e-5 5e-5])
set(gca, 'ydir', 'reverse', 'XAxisLocation', 'bottom', 'tickdir', 'out')
hold on
[c,h]=contour(Argo.latitude,Argo.pressure,SA_mean_E',[34.2,34.7,34.8,35.8],'k');
clabel(c,h);
[c,h]=contour(Argo.latitude,Argo.pressure,GM_mean_E',[25,25.8,26,26.2,26.8,27:0.1:28],'w');
clabel(c,h,'color','w');
ylabel('Pressure (dbar)');
xlabel('Latitude');
%set(gca, 'YTickLabel', {});
text(0.02, 1.12, '(g) $T''|_{\rm heave}$', 'Units', 'normalized', ...
     'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', ...
     'FontSize', 11, 'FontWeight', 'bold', 'Interpreter', 'latex');

% (h) Salinity pattern (Heave)
subplot('Position',[0.3,0.07,0.2,0.26]);
pcolor(lat(i_lat), Argo.pressure(i_d), EOF_maps_reshaped(:,1:58,1)'.*std(PCs(:,1)));
shading interp;
colorbar;
cmap_full = cmocean('balance', 36);
cmap_full = cmap_full([1:17, 20:36], :);   
colormap(gca, cmap_full);
caxis([-5e-5 5e-5])
set(gca, 'ydir', 'reverse', 'XAxisLocation', 'bottom', 'tickdir', 'out')
hold on
[c,h]=contour(Argo.latitude,Argo.pressure,SA_mean_E',[34.2,34.7,34.8,35.8],'k');
clabel(c,h);
[c,h]=contour(Argo.latitude,Argo.pressure,GM_mean_E',[25,25.8,26,26.2,26.8,27:0.1:28],'w');
clabel(c,h,'color','w');
ylabel('');
xlabel('Latitude');
set(gca, 'YTickLabel', {});
text(0.02, 1.12, '(h) $S''|_{\rm heave}$', 'Units', 'normalized', ...
     'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', ...
     'FontSize', 11, 'FontWeight', 'bold', 'Interpreter', 'latex');

% (i) PC time series (Heave)
subplot('Position',[0.58,0.1,0.4,0.18]);
plot(time, PCs(:,1)./std(PCs(:,1)), [pc1_color '-'], 'LineWidth', 1.5);
hold on
text(0.02, 1.12, '(i) $PC1|_{\rm heave}$', 'Units', 'normalized', ...
     'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', ...
     'FontSize', 11, 'FontWeight', 'bold', 'Interpreter', 'latex');
subtitle(['Explained variance: ',num2str(round(latent(1)/sum(latent)*100,1)),'%'])
xlabel('Time');
ylabel('Amplitude');
ylim([-2.5 2.5])
set(gca, 'tickdir', 'out')

pc1=PCs(:,1)./std(PCs(:,1));
pc1=(pc1-nanmean(pc1))./std(pc1);

dmi_smooth = movmean((DMI.Index_Argo-nanmean(DMI.Index_Argo))./std(DMI.Index_Argo,'omitnan'),12,1);
plot(DMI.Dates_Argo, dmi_smooth, 'Color', dmi_color, 'LineWidth', 1.5);
correlation_dmi = corr(pc1, dmi_smooth, 'rows', 'complete');
text(DMI.Dates_Argo(10), -2.2, ['Corr with DMI: ', num2str(round(correlation_dmi,2))]);

enso_smooth = movmean((nino34.Index_Argo-nanmean(nino34.Index_Argo))./std(nino34.Index_Argo,'omitnan'),12,1);
plot(nino34.Dates_Argo, enso_smooth, [enso_color '--'], 'LineWidth', 1.5);
correlation_enso = corr(pc1, enso_smooth, 'rows', 'complete');
text(nino34.Dates_Argo(120), -2.2, ['Corr with ENSO: ', num2str(round(correlation_enso,2))]);

plot(time, zeros(size(time)), 'Color', [0.5 0.5 0.5], 'LineStyle', '--', 'LineWidth', 0.8);
legend('PC1','DMI','ENSO', 'Location', 'best');
hold off

%% figuer 6 WSIO EOF
%% ========== Colour definitions (your own) ==========
pc1_color = 'k';            % PC1: black solid
dmi_color = [0.85 0.6 0.1]; % DMI: yellow-orange (solid)
enso_color = 'b';           % ENSO: blue dashed

% WSIO no smoothed EOF (40-70°E)
lo_range=[40 70];           % <<< CHANGED from [90 120] to [40 70]
lat=Argo.latitude;
la_range=[-65 0];
i_lat=find(Argo.latitude >= -65 & Argo.latitude <=0);
i_lo=find(Argo.longitude >=lo_range(1) & Argo.longitude <=lo_range(2));
i_d=find(Argo.pressure >= 0);
beta_sa=beta.*Argo.SA_month_all;
sa_mean=squeeze(nanmean(beta_sa(i_lo,i_lat,i_d,:),1));
sa_mean=movmean(deseason(sa_mean,Argo.t2),12,3);
alpha_ct=alpha.*Argo.CT_month_all;
ct_mean=squeeze(nanmean(alpha_ct(i_lo,i_lat,i_d,:),1));
ct_mean=movmean(deseason(ct_mean,Argo.t2),12,3);

[ny, nz, nt] = size(sa_mean);
sa_mean_reshape = reshape(sa_mean, ny*nz, nt);
ct_mean_reshape = reshape(ct_mean, ny*nz, nt); 
ct_sa_reshape=[sa_mean_reshape; ct_mean_reshape];
valid = all(~isnan(ct_sa_reshape), 2);
ct_sa_reshape_valid = ct_sa_reshape(valid, :);

ct_sa_anom = detrend(ct_sa_reshape_valid', 'linear');
[EOFs, PCs, latent] = pca(ct_sa_anom);

EOF_maps = nan(ny*nz, size(EOFs, 2));
EOF_maps(valid, :) = EOFs;
EOF_maps_reshaped = reshape(EOF_maps, ny, nz*2, []);

% plot total eof on depth level T|z S|z
figure('Position',[300,50,800,700])

% ========== Row 1: Total on Pressure Levels (a-c) ==========
% (a) Temperature spatial pattern
subplot('Position',[0.07,0.7,0.2,0.26]);
pcolor(lat(i_lat), Argo.pressure(i_d), EOF_maps_reshaped(:,59:end,1)'.*std(PCs(:,1)));
shading interp;
cmap_full = cmocean('balance', 36);
cmap_full = cmap_full([1:17, 20:36], :);   
colormap(gca, cmap_full);
caxis([-5e-5 5e-5])
set(gca, 'ydir', 'reverse', 'XAxisLocation', 'bottom', 'tickdir', 'out')
hold on
[c,h]=contour(Argo.latitude,Argo.pressure,SA_mean_E',[34.2,34.7,34.8,35.8],'k');
clabel(c,h);
[c,h]=contour(Argo.latitude,Argo.pressure,GM_mean_E',[25,25.8,26,26.2,26.8,27:0.1:28],'w');
clabel(c,h,'color','w');
ylabel('Pressure (dbar)');
xlabel('');
set(gca, 'XTickLabel', {});
% Label at top-left using text (no warning, perfect placement)
text(0.02, 1.13, '(a) $T''|_p$', 'Units', 'normalized', ...
     'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', ...
     'FontSize', 11, 'FontWeight', 'bold', 'Interpreter', 'latex');

% (b) Salinity spatial pattern
subplot('Position',[0.3,0.7,0.2,0.26]);
pcolor(lat(i_lat), Argo.pressure(i_d), EOF_maps_reshaped(:,1:58,1)'.*std(PCs(:,1)));
shading interp;
colorbar;
cmap_full = cmocean('balance', 36);
cmap_full = cmap_full([1:17, 20:36], :);   
colormap(gca, cmap_full);
caxis([-5e-5 5e-5])
set(gca, 'ydir', 'reverse', 'XAxisLocation', 'bottom', 'tickdir', 'out')
hold on
[c,h]=contour(Argo.latitude,Argo.pressure,SA_mean_E',[34.2,34.7,34.8,35.8],'k');
clabel(c,h);
[c,h]=contour(Argo.latitude,Argo.pressure,GM_mean_E',[25,25.8,26,26.2,26.8,27:0.1:28],'w');
clabel(c,h,'color','w');
ylabel('');
xlabel('');
set(gca, 'XTickLabel', {});
set(gca, 'YTickLabel', {});
text(0.02, 1.13, '(b) $S''|_p$', 'Units', 'normalized', ...
     'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', ...
     'FontSize', 11, 'FontWeight', 'bold', 'Interpreter', 'latex');

% (c) PC time series
subplot('Position',[0.58,0.75,0.4,0.18]);
plot(time, PCs(:,1)./std(PCs(:,1)), [pc1_color '-'], 'LineWidth', 1.5);
hold on
text(0.02, 1.13, '(c) $PC1|_p$', 'Units', 'normalized', ...
     'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', ...
     'FontSize', 11, 'FontWeight', 'bold', 'Interpreter', 'latex');
subtitle(['Explained variance: ',num2str(round(latent(1)/sum(latent)*100,1)),'%'])
xlabel('Time');
ylabel('Amplitude');
ylim([-2.5 2.5])
set(gca, 'tickdir', 'out')

pc1=PCs(:,1)./std(PCs(:,1));
pc1=(pc1-nanmean(pc1))./std(pc1);

% DMI: your own colour (yellow-orange)
dmi_smooth = movmean((DMI.Index_Argo-nanmean(DMI.Index_Argo))./std(DMI.Index_Argo,'omitnan'),12,1);
plot(DMI.Dates_Argo, dmi_smooth, 'Color', dmi_color, 'LineWidth', 1.5);
correlation_dmi = corr(pc1, dmi_smooth, 'rows', 'complete');
text(DMI.Dates_Argo(10), -2.2, ['Corr with DMI: ', num2str(round(correlation_dmi,2))]);

% ENSO: blue dashed
enso_smooth = movmean((nino34.Index_Argo-nanmean(nino34.Index_Argo))./std(nino34.Index_Argo,'omitnan'),12,1);
plot(nino34.Dates_Argo, enso_smooth, [enso_color '--'], 'LineWidth', 1.5);
correlation_enso = corr(pc1, enso_smooth, 'rows', 'complete');
text(nino34.Dates_Argo(120), -2.2, ['Corr with ENSO: ', num2str(round(correlation_enso,2))]);

% Grey zero line
plot(time, zeros(size(time)), 'Color', [0.5 0.5 0.5], 'LineStyle', '--', 'LineWidth', 0.8);
legend('PC1','DMI','ENSO', 'Location', 'best');
hold off


% ========== Row 2: Spice EOF (d-f) ==========
lo_range=[40 70];           % <<< CHANGED from [90 120] to [40 70]
lat=Argo.latitude;
la_range=[-65 0];
i_lat=find(Argo.latitude >= la_range(1) & Argo.latitude <=la_range(2));
i_lo=find(Argo.longitude >=lo_range(1) & Argo.longitude <=lo_range(2));
beta_sa=beta.*Argo.s_spice_z;
sa_mean=squeeze(nanmean(beta_sa(i_lo,i_lat,:,:),1));
sa_mean=movmean(deseason(sa_mean,Argo.t2),12,3);
alpha_ct=alpha.*Argo.t_spice_z;
ct_mean=squeeze(nanmean(alpha_ct(i_lo,i_lat,:,:),1));
ct_mean=movmean(deseason(ct_mean,Argo.t2),12,3);

[ny, nz, nt] = size(sa_mean);
sa_mean_reshape = reshape(sa_mean, ny*nz, nt);
valid_sa = all(~isnan(sa_mean_reshape), 2);
sa_reshape_valid=sa_mean_reshape(valid_sa,:);
ct_mean_reshape = reshape(ct_mean, ny*nz, nt);
valid_ct = all(~isnan(ct_mean_reshape), 2);
ct_reshape_valid=ct_mean_reshape(valid_ct,:);
ct_sa_reshape=[sa_reshape_valid; ct_reshape_valid];

ct_sa_anom = detrend(ct_sa_reshape', 'linear');
[EOFs, PCs, latent] = pca(ct_sa_anom);

EOF_maps = nan(ny*nz*2, size(EOFs, 2));
EOF_maps([valid_ct,valid_sa], :) = EOFs;
EOF_maps_reshaped = reshape(EOF_maps, ny, nz*2, []);

% (d) Temperature pattern (Spice)
subplot('Position',[0.07,0.38,0.2,0.26]);
pcolor(lat(i_lat), Argo.pressure(i_d), -EOF_maps_reshaped(:,59:end,1)'.*std(PCs(:,1)));
shading interp;
cmap_full = cmocean('balance', 36);
cmap_full = cmap_full([1:17, 20:36], :);   
colormap(gca, cmap_full);
caxis([-5e-5 5e-5])
set(gca, 'ydir', 'reverse', 'XAxisLocation', 'bottom', 'tickdir', 'out')
hold on
[c,h]=contour(Argo.latitude,Argo.pressure,SA_mean_E',[34.2,34.7,34.8,35.8],'k');
clabel(c,h);
[c,h]=contour(Argo.latitude,Argo.pressure,GM_mean_E',[25,25.8,26,26.2,26.8,27:0.1:28],'w');
clabel(c,h,'color','w');
ylabel('Pressure (dbar)');
xlabel('');
set(gca, 'XTickLabel', {});
text(0.02, 1.13, '(d) $T''|_{\rm spice}$', 'Units', 'normalized', ...
     'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', ...
     'FontSize', 11, 'FontWeight', 'bold', 'Interpreter', 'latex');

% (e) Salinity pattern (Spice)
subplot('Position',[0.3,0.38,0.2,0.26]);
pcolor(lat(i_lat), Argo.pressure(i_d), -EOF_maps_reshaped(:,1:58,1)'.*std(PCs(:,1)));
shading interp;
colorbar;
cmap_full = cmocean('balance', 36);
cmap_full = cmap_full([1:17, 20:36], :);   
colormap(gca, cmap_full);
caxis([-5e-5 5e-5])
set(gca, 'ydir', 'reverse', 'XAxisLocation', 'bottom', 'tickdir', 'out')
hold on
[c,h]=contour(Argo.latitude,Argo.pressure,SA_mean_E',[34.2,34.7,34.8,35.8],'k');
clabel(c,h);
[c,h]=contour(Argo.latitude,Argo.pressure,GM_mean_E',[25,25.8,26,26.2,26.8,27:0.1:28],'w');
clabel(c,h,'color','w');
ylabel('');
xlabel('');
set(gca, 'XTickLabel', {});
set(gca, 'YTickLabel', {});
text(0.02, 1.13, '(e) $S''|_{\rm spice}$', 'Units', 'normalized', ...
     'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', ...
     'FontSize', 11, 'FontWeight', 'bold', 'Interpreter', 'latex');

% (f) PC time series (Spice)
subplot('Position',[0.58,0.4,0.4,0.18]);
plot(time, -PCs(:,1)./std(PCs(:,1)), [pc1_color '-'], 'LineWidth', 1.5);
hold on
text(0.02, 1.13, '(f) $PC1|_{\rm spice}$', 'Units', 'normalized', ...
     'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', ...
     'FontSize', 11, 'FontWeight', 'bold', 'Interpreter', 'latex');
subtitle(['Explained variance: ',num2str(round(latent(1)/sum(latent)*100,1)),'%'])
xlabel('Time');
ylabel('Amplitude');
ylim([-2.5 2.5])
set(gca, 'tickdir', 'out')

pc1=PCs(:,1)./std(PCs(:,1));
pc1=(pc1-nanmean(pc1))./std(pc1);

dmi_smooth = movmean((DMI.Index_Argo-nanmean(DMI.Index_Argo))./std(DMI.Index_Argo,'omitnan'),12,1);
plot(DMI.Dates_Argo, dmi_smooth, 'Color', dmi_color, 'LineWidth', 1.5);
correlation_dmi = corr(pc1, dmi_smooth, 'rows', 'complete');
text(DMI.Dates_Argo(10), -2.2, ['Corr with DMI: ', num2str(round(correlation_dmi,2))]);

enso_smooth = movmean((nino34.Index_Argo-nanmean(nino34.Index_Argo))./std(nino34.Index_Argo,'omitnan'),12,1);
plot(nino34.Dates_Argo, enso_smooth, [enso_color '--'], 'LineWidth', 1.5);
correlation_enso = corr(pc1, enso_smooth, 'rows', 'complete');
text(nino34.Dates_Argo(120), -2.2, ['Corr with ENSO: ', num2str(round(correlation_enso,2))]);

plot(time, zeros(size(time)), 'Color', [0.5 0.5 0.5], 'LineStyle', '--', 'LineWidth', 0.8);
legend('PC1','DMI','ENSO', 'Location', 'best');
hold off


% ========== Row 3: Heave EOF (g-i) ==========
lo_range=[40 70];           % <<< CHANGED from [90 120] to [40 70]
lat=Argo.latitude;
la_range=[-65 0];
i_lat=find(Argo.latitude >= la_range(1) & Argo.latitude <=la_range(2));
i_lo=find(Argo.longitude >=lo_range(1) & Argo.longitude <=lo_range(2));
beta_sa=beta.*Argo.s_heave_z;
sa_mean=squeeze(nanmean(beta_sa(i_lo,i_lat,:,:),1));
sa_mean=movmean(deseason(sa_mean,Argo.t2),12,3);
alpha_ct=alpha.*Argo.t_heave_z;
ct_mean=squeeze(nanmean(alpha_ct(i_lo,i_lat,:,:),1));
ct_mean=movmean(deseason(ct_mean,Argo.t2),12,3);

[ny, nz, nt] = size(sa_mean);
sa_mean_reshape = reshape(sa_mean, ny*nz, nt);
valid_sa = all(~isnan(sa_mean_reshape), 2);
sa_reshape_valid=sa_mean_reshape(valid_sa,:);
ct_mean_reshape = reshape(ct_mean, ny*nz, nt);
valid_ct = all(~isnan(ct_mean_reshape), 2);
ct_reshape_valid=ct_mean_reshape(valid_ct,:);
ct_sa_reshape=[sa_reshape_valid; ct_reshape_valid];

ct_sa_anom = detrend(ct_sa_reshape', 'linear');
[EOFs, PCs, latent] = pca(ct_sa_anom);

EOF_maps = nan(ny*nz*2, size(EOFs, 2));
EOF_maps([valid_ct,valid_sa], :) = EOFs;
EOF_maps_reshaped = reshape(EOF_maps, ny, nz*2, []);

% (g) Temperature pattern (Heave)
subplot('Position',[0.07,0.07,0.2,0.26]);
pcolor(lat(i_lat), Argo.pressure(i_d), -EOF_maps_reshaped(:,59:end,1)'.*std(PCs(:,1)));
shading interp;
cmap_full = cmocean('balance', 36);
cmap_full = cmap_full([1:17, 20:36], :);   
colormap(gca, cmap_full);
caxis([-5e-5 5e-5])
set(gca, 'ydir', 'reverse', 'XAxisLocation', 'bottom', 'tickdir', 'out')
hold on
[c,h]=contour(Argo.latitude,Argo.pressure,SA_mean_E',[34.2,34.7,34.8,35.8],'k');
clabel(c,h);
[c,h]=contour(Argo.latitude,Argo.pressure,GM_mean_E',[25,25.8,26,26.2,26.8,27:0.1:28],'w');
clabel(c,h,'color','w');
ylabel('Pressure (dbar)');
xlabel('Latitude');
text(0.02, 1.12, '(g) $T''|_{\rm heave}$', 'Units', 'normalized', ...
     'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', ...
     'FontSize', 11, 'FontWeight', 'bold', 'Interpreter', 'latex');

% (h) Salinity pattern (Heave)
subplot('Position',[0.3,0.07,0.2,0.26]);
pcolor(lat(i_lat), Argo.pressure(i_d), -EOF_maps_reshaped(:,1:58,1)'.*std(PCs(:,1)));
shading interp;
colorbar;
cmap_full = cmocean('balance', 36);
cmap_full = cmap_full([1:17, 20:36], :);   
colormap(gca, cmap_full);
caxis([-5e-5 5e-5])
set(gca, 'ydir', 'reverse', 'XAxisLocation', 'bottom', 'tickdir', 'out')
hold on
[c,h]=contour(Argo.latitude,Argo.pressure,SA_mean_E',[34.2,34.7,34.8,35.8],'k');
clabel(c,h);
[c,h]=contour(Argo.latitude,Argo.pressure,GM_mean_E',[25,25.8,26,26.2,26.8,27:0.1:28],'w');
clabel(c,h,'color','w');
ylabel('');
xlabel('Latitude');
set(gca, 'YTickLabel', {});
text(0.02, 1.12, '(h) $S''|_{\rm heave}$', 'Units', 'normalized', ...
     'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', ...
     'FontSize', 11, 'FontWeight', 'bold', 'Interpreter', 'latex');

% (i) PC time series (Heave)
subplot('Position',[0.58,0.1,0.4,0.18]);
plot(time, -PCs(:,1)./std(PCs(:,1)), [pc1_color '-'], 'LineWidth', 1.5);
hold on
text(0.02, 1.12, '(i) $PC1|_{\rm heave}$', 'Units', 'normalized', ...
     'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', ...
     'FontSize', 11, 'FontWeight', 'bold', 'Interpreter', 'latex');
subtitle(['Explained variance: ',num2str(round(latent(1)/sum(latent)*100,1)),'%'])
xlabel('Time');
ylabel('Amplitude');
ylim([-2.5 2.5])
set(gca, 'tickdir', 'out')

pc1=PCs(:,1)./std(PCs(:,1));
pc1=(pc1-nanmean(pc1))./std(pc1);

dmi_smooth = movmean((DMI.Index_Argo-nanmean(DMI.Index_Argo))./std(DMI.Index_Argo,'omitnan'),12,1);
plot(DMI.Dates_Argo, dmi_smooth, 'Color', dmi_color, 'LineWidth', 1.5);
correlation_dmi = corr(pc1, dmi_smooth, 'rows', 'complete');
text(DMI.Dates_Argo(10), -2.2, ['Corr with DMI: ', num2str(round(correlation_dmi,2))]);

enso_smooth = movmean((nino34.Index_Argo-nanmean(nino34.Index_Argo))./std(nino34.Index_Argo,'omitnan'),12,1);
plot(nino34.Dates_Argo, enso_smooth, [enso_color '--'], 'LineWidth', 1.5);
correlation_enso = corr(pc1, enso_smooth, 'rows', 'complete');
text(nino34.Dates_Argo(120), -2.2, ['Corr with ENSO: ', num2str(round(correlation_enso,2))]);

plot(time, zeros(size(time)), 'Color', [0.5 0.5 0.5], 'LineStyle', '--', 'LineWidth', 0.8);
legend('PC1','DMI','ENSO', 'Location', 'best');
hold off


%% Figure 7 time sereis of ESIO
%prepare the data
%% calculate regional time series (20-30S 300-500m) in ESIO
i_lo=find(Argo.longitude>=90 & Argo.longitude <=120);
i_la=find(Argo.latitude >= -30  & Argo.latitude <=-20);
i_d=find(Argo.pressure >=300 & Argo.pressure <=500);
Argo.CT_ano=Argo.CT_month_all-nanmean(Argo.CT_month_all,4);
Argo.SA_ano=Argo.SA_month_all-nanmean(Argo.SA_month_all,4);
alpha_r=squeeze(nanmean(nanmean(nanmean(alpha_mean(i_lo,i_la,i_d,:),1),2),3));
beta_r=squeeze(nanmean(nanmean(nanmean(beta_mean(i_lo,i_la,i_d,:),1),2),3));
CT_z=alpha_r*squeeze(nanmean(nanmean(nanmean(Argo.CT_ano(i_lo,i_la,i_d,:),1),2),3));
SA_z=beta_r*squeeze(nanmean(nanmean(nanmean(Argo.SA_ano(i_lo,i_la,i_d,:),1),2),3));
CT_spice=alpha_r*squeeze(nanmean(nanmean(nanmean(Argo.t_spice_z(i_lo,i_la,i_d,:),1),2),3));
SA_spice=beta_r*squeeze(nanmean(nanmean(nanmean(Argo.s_spice_z(i_lo,i_la,i_d,:),1),2),3));
CT_heave=alpha_r*squeeze(nanmean(nanmean(nanmean(Argo.t_heave_z(i_lo,i_la,i_d,:),1),2),3));
SA_heave=beta_r*squeeze(nanmean(nanmean(nanmean(Argo.s_heave_z(i_lo,i_la,i_d,:),1),2),3));

%set colorbar
c_z     = [0, 0, 0];            % Black
c_spice = [49, 130, 189]/255;   % Blue (spice)
c_heave = [230, 85, 13]/255;    % Orange (heave)
c_grey =[0.5, 0.5, 0.5];


%
%% Figure setup
figure('Position',[100, 100, 1200, 900])
tiledlayout(3,2, 'Padding','compact', 'TileSpacing','compact')

% Line styles for climate indices
dmi_style  = {'--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.5}; % Gray dashed
enso_style = {':',  'Color', [0.2 0.4 0.9], 'LineWidth', 1.5}; % Blue dotted

% Terms & colors
terms  = {'p','spice','heave'};
colors = {c_z, c_spice, c_heave};
titles = {'T','S'};   % Variable names: T and S

% Variables to loop
vars = {CT_z, SA_z, CT_spice, SA_spice, -CT_heave, -SA_heave};

% Panel labels (a) to (f) in row-major order
panel_labels = {'(a)', '(b)', '(c)', '(d)', '(e)', '(f)'};

for j = 1:2   % Loop through variables: CT (1), SA (2)
    for i = 1:3 % Loop through terms: p, spice, heave
        idx = (i-1)*2 + j;      % Index for vars array
        data = vars{idx};       % Extract data
        color = colors{i};
        prefix = titles{j};

        % ---- Normalise ----
        norm_data = (data - nanmean(data)) ./ std(data,'omitnan');
        ENSO_norm = (nino34.Index_Argo - nanmean(nino34.Index_Argo)) ./ std(nino34.Index_Argo,'omitnan');
        DMI_norm  = (DMI.Index_Argo - nanmean(DMI.Index_Argo)) ./ std(DMI.Index_Argo,'omitnan');

        % ---- Plot ----
        nexttile(idx)
        hold on
        plot(Argo.t2, movmean(norm_data,12,1), ':', 'Color', c_heave, 'LineWidth', 2)
        
        area(Argo.t2, norm_data, 'FaceColor', c_heave, 'FaceAlpha', 0.4, 'EdgeAlpha', 0);
        plot(Argo.t2, -movmean(ENSO_norm,12,'omitnan'), 'color', c_spice, 'LineWidth', 1.5)
        plot(Argo.t2, -movmean(DMI_norm,12,'omitnan'), 'color', c_grey, 'LineWidth', 1.5)

        % ---- Only calculate lag & regression for z & heave ----
        if any(strcmp(terms{i}, {'p','heave'}))
            % ---- Lag calculation ----
            maxLag = 12; % months
            [acor_enso, lag_enso] = xcorr(norm_data, ENSO_norm, maxLag, 'coeff');
            [acor_dmi,  lag_dmi]  = xcorr(norm_data, DMI_norm, maxLag, 'coeff');
            [~, Ienso] = max(abs(acor_enso));
            [~, Idmi]  = max(abs(acor_dmi));
            bestLag_ENSO = lag_enso(Ienso);
            bestLag_DMI  = lag_dmi(Idmi);

            % ---- Shift ENSO & DMI ----
            ENSO_lagged = [nan(bestLag_ENSO,1); ENSO_norm(1:end-bestLag_ENSO)];
            DMI_lagged  = [nan(bestLag_DMI,1); DMI_norm(1:end-bestLag_DMI)];

            % ---- Align valid points for regression ----
            valid = ~isnan(norm_data) & ~isnan(ENSO_lagged) & ~isnan(DMI_lagged);
            X = [ones(sum(valid),1), ENSO_lagged(valid), DMI_lagged(valid)];
            y = norm_data(valid);

            % ---- Regression ----
            [b, ~, resid] = regress(y, X);
            norm_removed = nan(size(norm_data));
            norm_removed(valid) = resid; % Residuals: ENSO & DMI removed

            % ---- Plot residual line ----
            plot(Argo.t2, movmean(norm_removed,12,1), 'k', 'LineWidth', 1.8)
            polyplot(Argo.t2, norm_removed, 'color', 'k', 'linewidth', 1)
        end

        polyplot(Argo.t2, norm_data, 'color', c_heave) % trend line for data before mode removed

        % ---- Correlations ----
        [r_dmi, p_dmi]   = corr(movmean(norm_data,12,1), movmean(DMI_norm,12,'omitnan'));
        [r_enso, p_enso] = corr(movmean(norm_data,12,1), movmean(ENSO_norm,12,'omitnan'));
        fprintf('ENSO coefficient = %.3f\n', r_enso)
        fprintf('DMI coefficient = %.3f\n',  r_dmi)

        % Multiple linear regression
        y  = movmean(norm_data,12,1);
        x1 = movmean(ENSO_norm,12,'omitnan');
        x2 = movmean(DMI_norm,12,'omitnan');
        X = [x1 x2];
        [b, bint, r, rint, stats] = regress(y, [ones(size(X,1),1) X]);
        R2 = stats(1);
        p_combined = stats(3);
        fprintf('Combined R² = %.3f, p = %.3f\n', R2, p_combined)
        fprintf('ENSO coefficient = %.3f\n', b(2))
        fprintf('DMI coefficient = %.3f\n',  b(3))

        text(Argo.t2(10), -2.5, sprintf('r(DMI) = %.2f, p = %.3f', r_dmi, p_dmi), 'FontSize', 9)
        text(Argo.t2(10), -3.0, sprintf('r(ENSO) = %.2f, p = %.3f', r_enso, p_enso), 'FontSize', 9)
        ylim([-4 4])
        ylabel('Standardized')
        if i == 3
            xlabel('Year')
        end

        % ---- TITLE: (a) $T'$ | p  etc. placed at top-left inside axes ----
        title_str = sprintf('%s $%s''$ | %s', panel_labels{idx}, prefix, terms{i});
        title(title_str, ...
              'Position', [0.02, 4.8], ...          % top-left inside axes
              'Units', 'normalized', ...
              'HorizontalAlignment', 'left', ...
              'VerticalAlignment', 'top', ...
              'Interpreter', 'latex', ...
              'FontSize', 11, 'FontWeight', 'bold')

        % ---- Legend ----
        if strcmp(terms{i}, 'spice')
            legend([titles{j},'|', terms{i},' 12 month smoothed'], ...
                   [titles{j},'|', terms{i}], '-ENSO', '-DMI', ...
                   'Location', 'northeast', 'Box', 'off')
        else
            legend([titles{j},'|', terms{i},' 12 month smoothed'], ...
                   [titles{j},'|', terms{i}], '-ENSO', '-DMI', ...
                   [titles{j},'|', terms{i},' removed'], ...
                   'Location', 'northeast', 'Box', 'off')
        end
        set(gca, 'FontSize', 10, 'TickDir', 'out')
    end
end


%% Figure 8 time sereis of WSIO
%prepare data
%% calculate regional time series (20-30S 300-500m) in ESIO
i_lo=find(Argo.longitude>=40 & Argo.longitude <=70);
i_la=find(Argo.latitude >= -30  & Argo.latitude <=-20);
i_d=find(Argo.pressure >=300 & Argo.pressure <=500);
Argo.CT_ano=Argo.CT_month_all-nanmean(Argo.CT_month_all,4);
Argo.SA_ano=Argo.SA_month_all-nanmean(Argo.SA_month_all,4);
alpha_r=squeeze(nanmean(nanmean(nanmean(alpha_mean(i_lo,i_la,i_d,:),1),2),3));
beta_r=squeeze(nanmean(nanmean(nanmean(beta_mean(i_lo,i_la,i_d,:),1),2),3));
CT_z=alpha_r*squeeze(nanmean(nanmean(nanmean(Argo.CT_ano(i_lo,i_la,i_d,:),1),2),3));
SA_z=beta_r*squeeze(nanmean(nanmean(nanmean(Argo.SA_ano(i_lo,i_la,i_d,:),1),2),3));
CT_spice=alpha_r*squeeze(nanmean(nanmean(nanmean(Argo.t_spice_z(i_lo,i_la,i_d,:),1),2),3));
SA_spice=beta_r*squeeze(nanmean(nanmean(nanmean(Argo.s_spice_z(i_lo,i_la,i_d,:),1),2),3));
CT_heave=alpha_r*squeeze(nanmean(nanmean(nanmean(Argo.t_heave_z(i_lo,i_la,i_d,:),1),2),3));
SA_heave=beta_r*squeeze(nanmean(nanmean(nanmean(Argo.s_heave_z(i_lo,i_la,i_d,:),1),2),3));

%set colorbar
c_z     = [0, 0, 0];            % Black
c_spice = [49, 130, 189]/255;   % Blue (spice)
c_heave = [230, 85, 13]/255;    % Orange (heave)
c_grey =[0.5, 0.5, 0.5];

%%
%% Figure setup
figure('Position',[100, 100, 1200, 900])
tiledlayout(3,2, 'Padding','compact', 'TileSpacing','compact')

% Line styles for climate indices
dmi_style  = {'--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.5}; % Gray dashed
enso_style = {':',  'Color', [0.2 0.4 0.9], 'LineWidth', 1.5}; % Blue dotted

% Terms & colors
terms  = {'p','spice','heave'};
colors = {c_z, c_spice, c_heave};
titles = {'T','S'};   % Variable names: T and S

% Variables to loop
vars = {CT_z, SA_z, CT_spice, SA_spice, -CT_heave, -SA_heave};

% Panel labels (a) to (f) in row-major order
panel_labels = {'(a)', '(b)', '(c)', '(d)', '(e)', '(f)'};

for j = 1:2   % Loop through variables: CT (1), SA (2)
    for i = 1:3 % Loop through terms: p, spice, heave
        idx = (i-1)*2 + j;      % Index for vars array
        data = vars{idx};       % Extract data
        color = colors{i};
        prefix = titles{j};

        % ---- Normalise ----
        norm_data = (data - nanmean(data)) ./ std(data,'omitnan');
        ENSO_norm = (nino34.Index_Argo - nanmean(nino34.Index_Argo)) ./ std(nino34.Index_Argo,'omitnan');
        DMI_norm  = (DMI.Index_Argo - nanmean(DMI.Index_Argo)) ./ std(DMI.Index_Argo,'omitnan');

        % ---- Plot ----
        nexttile(idx)
        hold on
        plot(Argo.t2, movmean(norm_data,12,1), ':', 'Color', c_heave, 'LineWidth', 2)
        
        area(Argo.t2, norm_data, 'FaceColor', c_heave, 'FaceAlpha', 0.4, 'EdgeAlpha', 0);
        plot(Argo.t2, -movmean(ENSO_norm,12,'omitnan'), 'color', c_spice, 'LineWidth', 1.5)
        plot(Argo.t2, -movmean(DMI_norm,12,'omitnan'), 'color', c_grey, 'LineWidth', 1.5)

        % ---- Only calculate lag & regression for z & heave ----
        if any(strcmp(terms{i}, {'p','heave'}))
            % ---- Lag calculation ----
            maxLag = 12; % months
            [acor_enso, lag_enso] = xcorr(norm_data, ENSO_norm, maxLag, 'coeff');
            [acor_dmi,  lag_dmi]  = xcorr(norm_data, DMI_norm, maxLag, 'coeff');
            [~, Ienso] = max(abs(acor_enso));
            [~, Idmi]  = max(abs(acor_dmi));
            bestLag_ENSO = lag_enso(Ienso);
            bestLag_DMI  = lag_dmi(Idmi);

            % ---- Shift ENSO & DMI ----
            ENSO_lagged = [nan(bestLag_ENSO,1); ENSO_norm(1:end-bestLag_ENSO)];
            DMI_lagged  = [nan(bestLag_DMI,1); DMI_norm(1:end-bestLag_DMI)];

            % ---- Align valid points for regression ----
            valid = ~isnan(norm_data) & ~isnan(ENSO_lagged) & ~isnan(DMI_lagged);
            X = [ones(sum(valid),1), ENSO_lagged(valid), DMI_lagged(valid)];
            y = norm_data(valid);

            % ---- Regression ----
            [b, ~, resid] = regress(y, X);
            norm_removed = nan(size(norm_data));
            norm_removed(valid) = resid; % Residuals: ENSO & DMI removed

            % ---- Plot residual line ----
            plot(Argo.t2, movmean(norm_removed,12,1), 'k', 'LineWidth', 1.8)
            polyplot(Argo.t2, norm_removed, 'color', 'k', 'linewidth', 1)
        end

        polyplot(Argo.t2, norm_data, 'color', c_heave) % trend line for data before mode removed

        % ---- Correlations ----
        [r_dmi, p_dmi]   = corr(movmean(norm_data,12,1), movmean(DMI_norm,12,'omitnan'));
        [r_enso, p_enso] = corr(movmean(norm_data,12,1), movmean(ENSO_norm,12,'omitnan'));
        fprintf('ENSO coefficient = %.3f\n', r_enso)
        fprintf('DMI coefficient = %.3f\n',  r_dmi)

        % Multiple linear regression
        y  = movmean(norm_data,12,1);
        x1 = movmean(ENSO_norm,12,'omitnan');
        x2 = movmean(DMI_norm,12,'omitnan');
        X = [x1 x2];
        [b, bint, r, rint, stats] = regress(y, [ones(size(X,1),1) X]);
        R2 = stats(1);
        p_combined = stats(3);
        fprintf('Combined R² = %.3f, p = %.3f\n', R2, p_combined)
        fprintf('ENSO coefficient = %.3f\n', b(2))
        fprintf('DMI coefficient = %.3f\n',  b(3))

        text(Argo.t2(10), -2.5, sprintf('r(DMI) = %.2f, p = %.3f', r_dmi, p_dmi), 'FontSize', 9)
        text(Argo.t2(10), -3.0, sprintf('r(ENSO) = %.2f, p = %.3f', r_enso, p_enso), 'FontSize', 9)
        ylim([-4 4])
        ylabel('Standardized')
        if i == 3
            xlabel('Year')
        end

        % ---- TITLE: (a) $T'$ | p  etc. placed at top-left inside axes ----
        title_str = sprintf('%s $%s''$ | %s', panel_labels{idx}, prefix, terms{i});
        title(title_str, ...
              'Position', [0.02, 4.8], ...          % top-left inside axes
              'Units', 'normalized', ...
              'HorizontalAlignment', 'left', ...
              'VerticalAlignment', 'top', ...
              'Interpreter', 'latex', ...
              'FontSize', 11, 'FontWeight', 'bold')

        % ---- Legend ----
        if strcmp(terms{i}, 'spice')
            legend([titles{j},'|', terms{i},' 12 month smoothed'], ...
                   [titles{j},'|', terms{i}], '-ENSO', '-DMI', ...
                   'Location', 'northeast', 'Box', 'off')
        else
            legend([titles{j},'|', terms{i},' 12 month smoothed'], ...
                   [titles{j},'|', terms{i}], '-ENSO', '-DMI', ...
                   [titles{j},'|', terms{i},' removed'], ...
                   'Location', 'northeast', 'Box', 'off')
        end
        set(gca, 'FontSize', 10, 'TickDir', 'out')
    end
end


%% Figure 9 T anomaly along 26.8
% --- Extract temperature on gamma=26.8 surface ---
land = shaperead('landareas.shp', 'UseGeoCoords', true);
load coastlines.mat
i_g = find(abs(Argo.gamma_grid - 26.8) < 1e-3, 1); % safer than == 
T_iso = squeeze(Argo.t_spice(:,:,i_g,:)); % [lon × lat × time]

% --- Time info ---
time = Argo.t2; % assume datetime format
years = year(time);
months = month(time);

% --- Calculate annual means ---
unique_years = unique(years);
nYears = numel(unique_years);

T_annual = nan([size(T_iso,1), size(T_iso,2), nYears]);
for i = 1:nYears
    idx = years == unique_years(i);
    T_annual(:,:,i) = mean(T_iso(:,:,idx), 3, 'omitnan');
end

% --- Define anomaly relative to 2014 ---
i_ref = find(unique_years == 2014);
T_ref = T_annual(:,:,i_ref); %2004
T_mean=squeeze(nanmean(T_annual,3));
T_anom = T_annual - T_mean; % anomaly to mean
%T_anom = T_annual - T_ref; % anomaly to 2004
%% figure plot
% --- Plot settings ---
nRow = 4; nCol = 5;
x_margin = 0.06; y_margin = 0.06;
x_spacing = 0.015; y_spacing = 0.05;
subplot_width  = (1 - 2*x_margin - (nCol-1)*x_spacing)/nCol;
subplot_height = (1 - 2*y_margin - (nRow-1)*y_spacing)/nRow;

letters = 'abcdefghijklmnopqrst'; % 20 panels (a–t)

figure('Position',[100 100 1400 900]);

plot_idx = 0; % counter for subplots (skip 2014)

for i = 1:21
    yr = 2003 + i; % years 2004–2024
    if yr == 2014
        continue; % skip this year
    end
    
    plot_idx = plot_idx + 1; % increment subplot index
    idx = find(unique_years == yr);

    row = ceil(plot_idx/nCol); % subplot row index
    col = mod(plot_idx-1,nCol) + 1;

    left   = x_margin + (col-1)*(subplot_width + x_spacing);
    bottom = 1 - y_margin - row*subplot_height - (row-1)*y_spacing;

    ax = axes('Position',[left bottom subplot_width subplot_height]); 
    hold on; box on;
    
    pcolor(Argo.longitude, Argo.latitude, T_anom(:,:,idx)');
    shading flat; axis tight;
    caxis([-0.5 0.5]); % symmetric anomaly range
    hold on
    geoshow(land, 'FaceColor', [0.85 0.85 0.85], 'EdgeColor', 'none');
        plot(coastlon, coastlat, 'k-', 'LineWidth', 0.6);
    title(sprintf('(%s) %d', letters(plot_idx), yr), 'FontWeight','normal');
    set(gca,'TickDir','out','Layer','top');
    ylim([-50 0])
    xlim([20 150]);
    % Hide ticks for internal plots
    if row < nRow
        set(gca,'XTickLabel',[]);
    else
        xlabel('Longitude');
    end
    if col > 1
        set(gca,'YTickLabel',[]);
    else
        ylabel('Latitude');
    end
end

% Shared colorbar
h = colorbar('Position',[0.92 0.1 0.02 0.8]); 
ylabel(h,'Temperature anomaly (°C)');

% Adjust colormap (remove center white gap if needed)
cmap = cmocean('balance',20);
cmap = cmap([1:9,12:20],:);
colormap(cmap);



%% Figure 10 hovmoller figure along 26.8
%% Figure setup
figure('Position',[100 100 1000 900]); 

nRow = 3; % now 3 rows: upper 2 for longitude ranges, lower 1 for 20-40S avg
nCol = 3;

x_margin = 0.07; y_margin = 0.06;
x_spacing = 0.05; y_spacing = 0.06;
subplot_width  = (1 - 2*x_margin - (nCol-1)*x_spacing)/nCol - 0.01;
subplot_height = (1 - 2*y_margin - (nRow-1)*y_spacing)/nRow-0.005;

letters = 'abcdefghi'; % 9 panels (a-i)
lon_ranges = {[40 70],[90 120]};
terms = {'t_z_n','t_spice','t_heave'};
titles_text = {'T''|p','T''|spice','T''|heave'};
scale_factor = [1 1 -1]; % negate heave
i_g = find(Argo.gamma_grid==26.8);

% --- Upper 2 rows: longitude bands ---
for r = 1:2
    i_lon = find(Argo.longitude >= lon_ranges{r}(1) & Argo.longitude <= lon_ranges{r}(2));
    
    for c = 1:nCol
        left   = x_margin + (c-1)*(subplot_width + x_spacing);
        bottom = 1 - y_margin - r*subplot_height - (r-1)*y_spacing;
        ax = axes('Position',[left bottom subplot_width subplot_height]); hold on; box on;
        
        % extract and smooth data
        var_data = squeeze(nanmean(eval(['Argo.' terms{c} '(i_lon,:,i_g,:)']),1));
        var_data = movmean(scale_factor(c)*var_data,12,2);
        
        pcolor(Argo.latitude, Argo.t2, var_data'); shading flat;
        caxis([-0.4 0.4]);
        cmap_full = cmocean('balance', 36);
        cmap_full = cmap_full([1:17, 20:36], :);   % 15 blue + 15 red, perfectly balanced

        colormap(gca, cmap_full);
        set(gca,'tickdir','out');
        xlim([-50 0])
        
        % Add °S to latitude ticks
        xt = get(gca, 'XTick');
        xt_labels = cell(size(xt));
        for i = 1:length(xt)
            if xt(i) == 0
                xt_labels{i} = '0°';
            else
                xt_labels{i} = sprintf('%d°S', abs(xt(i)));
            end
        end
        set(gca, 'XTickLabel', xt_labels);
        
        % Labels
        if r==2, xlabel('Latitude (^\circS)'); end
        if c==1, ylabel('Time'); end
        
        % Title with panel label
        title(sprintf('(%s) %s(%d-%d°E)', letters((r-1)*3+c), titles_text{c}, lon_ranges{r}(1), lon_ranges{r}(2)), 'FontWeight','normal');
    end
end

% --- Lower row: latitude-averaged 20-40S ---
i_lat_20_40 = find(Argo.latitude <= -20 & Argo.latitude >= -40);
r = 3; % third row

for c = 1:nCol
    left   = x_margin + (c-1)*(subplot_width + x_spacing);
    bottom = 1 - y_margin - r*subplot_height - (r-1)*y_spacing;
    ax = axes('Position',[left bottom-0.02 subplot_width subplot_height]); hold on; box on;
    
    % extract and smooth data averaged over latitude
    var_data = squeeze(nanmean(eval(['Argo.' terms{c} '(:,i_lat_20_40,i_g,:)']),2));
    var_data = movmean(scale_factor(c)*var_data,12,2);
    
    pcolor(Argo.longitude, Argo.t2, var_data'); shading flat;
    caxis([-0.4 0.4]);
    cmap_full = cmocean('balance', 36);
    cmap_full = cmap_full([1:17, 20:36], :);   % 15 blue + 15 red, perfectly balanced
    xlim([20 140])

    colormap(gca, cmap_full);
    set(gca,'tickdir','out');
    
    % Add °E to longitude ticks
    xt = get(gca, 'XTick');
    xt_labels = arrayfun(@(x) sprintf('%d°E', x), xt, 'UniformOutput', false);
    set(gca, 'XTickLabel', xt_labels);
    
    xlabel('Longitude (^\circE)');
    if c==1, ylabel('Time'); end
    
    % Title
    title(sprintf('(%s) %s(%d-%d°S)', letters(6+c), titles_text{c},20,40), 'FontWeight','normal');
end

% --- Shared colorbar ---
h = colorbar('Position',[0.92 0.08 0.02 0.84]);
h.TickDirection = 'out';  % <--- Make colorbar ticks point outward
ylabel(h,'Temperature anomaly (°C)');

cmap_full = cmocean('balance', 36);
cmap_full = cmap_full([1:17, 20:36], :);   % 15 blue + 15 red, perfectly balanced
colormap(gca, cmap_full);
caxis([-0.4 0.4]);


%% Figure 11 density trend and surface density change
  %% plot together
    figure('Position',[100 50 800 800]);

% Define positions: [left, bottom, width, height]
pos.a = [0.08 0.58 0.38 0.35];  % upper left
pos.b = [0.55 0.58 0.38 0.35];  % upper right
pos.c = [0.08 0.10 0.85 0.35];  % bottom (wide)

%% (a) Isopycnal trend 40–70E
axes('Position',pos.a);
i_lo = find(Argo.longitude >=40 & Argo.longitude <=70);

% --- CHANGED TO ESS TREND ---
% [iso_tr, iso_tr_p] = trend(squeeze(nanmean(Argo.gamma_n(i_lo,:,:,:),1)),12,'omitnan');
[iso_tr, iso_tr_p] = calc_trend_ess(squeeze(nanmean(Argo.gamma_n(i_lo,:,:,:),1)));

pcolor(Argo.latitude, Argo.pressure, iso_tr'); shading flat
hold on
% 注意：原代码 >0.05 表示标记“不显著”区域，若需标记显著区请改为 <0.05
StatisticallySignificant = iso_tr_p > 0.05;
stipple(Argo.latitude, Argo.pressure, StatisticallySignificant', ...
        'color', [0.4 0.4 0.4], 'density', 70, 'marker', '.', 'markersize', 4);

[c,h] = contour(Argo.latitude, Argo.pressure, squeeze(nanmean(Argo.gm(i_lo,:,:),1))', ...
                [25,25.8,26,26.2,26.8,27:0.1:28], 'k','LineWidth',1.2);
clabel(c,h,'color','k');

caxis([-0.006 0.006])
cmap_full = cmocean('balance', 36);
cmap_full = cmap_full([1:17, 20:36], :);
colormap(gca, cmap_full);

ylabel('Pressure (dbar)')
set(gca,'YDir','reverse', 'TickDir','out', 'FontSize',10);
% 为纬度轴（x轴）添加 °S 标签
lat_ticks = -50:10:0;
set(gca, 'XTick', lat_ticks, ...
    'XTickLabel', arrayfun(@(x) sprintf('%d°S', abs(x)), lat_ticks, 'UniformOutput', false));
xlabel('Latitude')
title('(a) Isopycnal trend western SIO')

%% (b) Isopycnal trend 90–120E
axes('Position',pos.b);
i_lo = find(Argo.longitude >=90 & Argo.longitude <=120);

% --- CHANGED TO ESS TREND ---
% [iso_tr, iso_tr_p] = trend(squeeze(nanmean(Argo.gamma_n(i_lo,:,:,:),1)),12,'omitnan');
[iso_tr, iso_tr_p] = calc_trend_ess(squeeze(nanmean(Argo.gamma_n(i_lo,:,:,:),1)));

pcolor(Argo.latitude, Argo.pressure, iso_tr'); shading flat
hold on
StatisticallySignificant = iso_tr_p > 0.05;
stipple(Argo.latitude, Argo.pressure, StatisticallySignificant', ...
        'color', [0.4 0.4 0.4], 'density', 70, 'marker', '.', 'markersize', 4);

[c,h] = contour(Argo.latitude, Argo.pressure, squeeze(nanmean(Argo.gm(i_lo,:,:),1))', ...
                [25,25.8,26,26.2,26.8,27:0.1:28], 'k','LineWidth',1.2);
clabel(c,h,'color','k');

caxis([-0.006 0.006])
cmap_full = cmocean('balance', 36);
cmap_full = cmap_full([1:17, 20:36], :);
colormap(gca, cmap_full);

cb = colorbar;
cb.Label.String = 'Neutral density trend (kg m^{-3} yr^{-1})';
cb.TickDirection = 'out'; % 确保 colorbar 刻度向外

set(gca,'YDir','reverse', 'YTickLabel', [], 'TickDir','out', 'FontSize',10);
% 纬度轴标签
set(gca, 'XTick', lat_ticks, ...
    'XTickLabel', arrayfun(@(x) sprintf('%d°S', abs(x)), lat_ticks, 'UniformOutput', false));
xlabel('Latitude')
title('(b) Isopycnal trend eastern SIO')

%% (c) Surface neutral density difference
axes('Position',pos.c);
surf_diff = squeeze(nanmean(Argo.gamma_n(:,:,1,193:252),4))' ...
          - squeeze(nanmean(Argo.gamma_n(:,:,1,1:60),4))';

p1 = pcolor(Argo.longitude, Argo.latitude, surf_diff); shading flat
hold on

% 2020–2024 等密度线（绿色虚线）
h1 = contour(Argo.longitude, Argo.latitude, ...
          squeeze(nanmean(Argo.gamma_n(:,:,1,193:252),4))', ...
          [21,22,23,24,25,25.8,26.2,26.8,27,27.2], 'g:','LineWidth',1.8);

geoshow(land, 'FaceColor', [0.85 0.85 0.85], 'EdgeColor', 'none');
plot(coastlon, coastlat, 'k-', 'LineWidth', 0.6);

% 2004–2008 等密度线（黑色实线）
h2 = contour(Argo.longitude, Argo.latitude, ...
          squeeze(nanmean(Argo.gamma_n(:,:,1,1:60),4))', ...
          [21,22,23,24,25,25.8,26.2,26.8,27,27.2], 'k-','LineWidth',1.4);
clabel(h2, 'color','k');

caxis([-0.3 0.3])
cmap_full = cmocean('balance', 36);
cmap_full = cmap_full([1:17, 20:36], :);
colormap(gca, cmap_full);

ch = colorbar;
ch.Label.String = 'Neutral density difference (kg m^{-3})';
ch.TickDirection = 'out';

% 使用 apply_geo_ticks 统一经纬度刻度
apply_geo_ticks(gca, 20:20:140, -50:10:0, true, true);
set(gca, 'TickDir','out', 'FontSize',10, 'Layer','top');
xlim([20 150]); ylim([-50 0]);

xlabel('Longitude')
ylabel('Latitude')
title('(c) Surface neutral density difference: 2020–2024 minus 2004–2008')

legend('Difference','2020-2024','2004-2008');