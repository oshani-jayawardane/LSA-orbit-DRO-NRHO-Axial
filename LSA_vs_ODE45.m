clc;
clear;
close all;

% ============================================================
% User selection
% Options: 'DRO', 'NRHO', 'Axial'
% ============================================================

Orbit = 'Axial';

ode_file = sprintf('datasets/ode45_Single_Orbit_%s_Dataset.mat', Orbit);
lsa_file = sprintf('LSA/LSA_%s_Table.mat', Orbit);

results_folder = 'results';

if ~exist(results_folder, 'dir')
    mkdir(results_folder);
end

% ============================================================
% Constants
% ============================================================

MU = 0.0121505856;
lstar = 384400;
tstar = 3.751892968837575e+05;

% ============================================================
% Load ODE45 reference dataset
% Dataset format: [t, x, y, z, vx, vy, vz, ax, ay, az]
% ============================================================

ref = load(ode_file);
Dataset = ref.Dataset;

if isfield(ref, 'tstar')
    tstar = ref.tstar;
end

t_ref = Dataset(:,1);
X_ref = Dataset(:,2:10);

pos_ref = X_ref(:,1:3);
vel_ref = X_ref(:,4:6);
acc_ref = X_ref(:,7:9);

% Time variables
% tau_lsa is used for polynomial evaluation
% tau_plot is used only for plotting from 0 to 1
tau_lsa = 2*(t_ref - t_ref(1)) / (t_ref(end) - t_ref(1)) - 1;
tau_plot = (t_ref - t_ref(1)) / (t_ref(end) - t_ref(1));

% ============================================================
% Load LSA coefficient table
% Final_Table has:
% Component | BestDegree | c0 | c1 | ... | c20
% ============================================================

lsa = load(lsa_file);
Final_Table = lsa.Final_Table;

components = string(Final_Table.Component);
best_degrees = Final_Table.BestDegree;

num_points = length(t_ref);
num_components = 9;

X_lsa = zeros(num_points, num_components);

expected_components = ["x","y","z","vx","vy","vz","ax","ay","az"];

% ============================================================
% Reconstruct LSA values from polynomial coefficients
% polynomial: c0 + c1*tau + c2*tau^2 + ...
% ============================================================

for j = 1:num_components

    comp_name = expected_components(j);

    row_idx = find(components == comp_name, 1);

    if isempty(row_idx)
        error('Component %s not found in Final_Table.', comp_name);
    end

    degree = best_degrees(row_idx);

    coeffs = table2array(Final_Table(row_idx, 3:end));
    coeffs = coeffs(1:degree+1);
    coeffs = coeffs(~isnan(coeffs));

    y_lsa = zeros(num_points,1);

    for k = 1:length(coeffs)
        y_lsa = y_lsa + coeffs(k) * tau_lsa.^(k-1);
    end

    X_lsa(:,j) = y_lsa;
end

pos_lsa = X_lsa(:,1:3);
vel_lsa = X_lsa(:,4:6);
acc_lsa = X_lsa(:,7:9);

% ============================================================
% Error computation
% Pointwise vector errors
% ============================================================

pos_error_nd = pos_lsa - pos_ref;
vel_error_nd = vel_lsa - vel_ref;
acc_error_nd = acc_lsa - acc_ref;

pos_error_km = sqrt(sum(pos_error_nd.^2, 2)) * lstar;
vel_error_kms = sqrt(sum(vel_error_nd.^2, 2)) * (lstar / tstar);
acc_error_kms2 = sqrt(sum(acc_error_nd.^2, 2)) * (lstar / tstar^2);

position_rmse_km = sqrt(mean(sum(pos_error_nd.^2, 2))) * lstar;
velocity_rmse_kms = sqrt(mean(sum(vel_error_nd.^2, 2))) * (lstar / tstar);
acceleration_rmse_kms2 = sqrt(mean(sum(acc_error_nd.^2, 2))) * (lstar / tstar^2);

fprintf('\n============================================================\n');
fprintf('LSA vs ODE45 Reference Error Summary\n');
fprintf('============================================================\n');

fprintf('Position error max  : %.6e km\n', max(pos_error_km));
fprintf('Position error mean : %.6e km\n', mean(pos_error_km));
fprintf('Position RMSE       : %.6e km\n', position_rmse_km);

fprintf('Velocity error max  : %.6e km/s\n', max(vel_error_kms));
fprintf('Velocity error mean : %.6e km/s\n', mean(vel_error_kms));
fprintf('Velocity RMSE       : %.6e km/s\n', velocity_rmse_kms);

fprintf('Acceleration error max  : %.6e km/s^2\n', max(acc_error_kms2));
fprintf('Acceleration error mean : %.6e km/s^2\n', mean(acc_error_kms2));
fprintf('Acceleration RMSE       : %.6e km/s^2\n', acceleration_rmse_kms2);

% ============================================================
% Convert trajectory to km for plotting
% ============================================================

x_ref_km = pos_ref(:,1) * lstar;
y_ref_km = pos_ref(:,2) * lstar;
z_ref_km = pos_ref(:,3) * lstar;

x_lsa_km = pos_lsa(:,1) * lstar;
y_lsa_km = pos_lsa(:,2) * lstar;
z_lsa_km = pos_lsa(:,3) * lstar;

all_x = [x_ref_km; x_lsa_km];
all_y = [y_ref_km; y_lsa_km];
all_z = [z_ref_km; z_lsa_km];

cx = (min(all_x) + max(all_x)) / 2;
cy = (min(all_y) + max(all_y)) / 2;
cz = (min(all_z) + max(all_z)) / 2;

span_x = max(all_x) - min(all_x);
span_y = max(all_y) - min(all_y);
span_z = max(all_z) - min(all_z);

pad = 1.10;
span_2d = max(span_x, span_y) * pad;
span_3d = max([span_x, span_y, span_z]) * pad;

min_span = 1e4;
span_2d = max(span_2d, min_span);
span_3d = max(span_3d, min_span);

% ============================================================
% 2D trajectory comparison
% ============================================================

figure;
DrawMoonCR3BP(MU, 1741, 1741)
hold on

h1 = plot3(x_ref_km, y_ref_km, z_ref_km, 'b', 'LineWidth', 1.5);
h2 = plot3(x_lsa_km, y_lsa_km, z_lsa_km, 'r', 'LineWidth', 1.5);

plot3(x_ref_km(1), y_ref_km(1), z_ref_km(1), 'ko', 'MarkerSize', 8, 'LineWidth', 1.5);
plot3(x_ref_km(end), y_ref_km(end), z_ref_km(end), 'kx', 'MarkerSize', 10, 'LineWidth', 1.5);

set(gcf,'color','w')
axis equal
grid on
view(2)

xlim(cx + [-span_2d, span_2d]/2)
ylim(cy + [-span_2d, span_2d]/2)

xlabel('X [km]')
ylabel('Y [km]')
zlabel('Z [km]')
title(sprintf('Trajectory Approximation (%s Orbit)', Orbit))

legend([h1 h2], {'ODE45', 'LSA'}, 'Location', 'best')

filename = fullfile(results_folder, sprintf('%s_2D_Trajectory_Comparison.png', Orbit));
exportgraphics(gcf, filename, 'Resolution', 300);

% ============================================================
% 3D trajectory comparison
% ============================================================

figure;
DrawMoonCR3BP(MU, 1741, 1741)
hold on

h1 = plot3(x_ref_km, y_ref_km, z_ref_km, 'b', 'LineWidth', 1.5);
h2 = plot3(x_lsa_km, y_lsa_km, z_lsa_km, 'r', 'LineWidth', 1.5);

plot3(x_ref_km(1), y_ref_km(1), z_ref_km(1), 'ko', 'MarkerSize', 8, 'LineWidth', 1.5);
plot3(x_ref_km(end), y_ref_km(end), z_ref_km(end), 'kx', 'MarkerSize', 10, 'LineWidth', 1.5);

set(gcf,'color','w')
axis equal
grid on
view(35,25)

xlim(cx + [-span_3d, span_3d]/2)
ylim(cy + [-span_3d, span_3d]/2)
zlim(cz + [-span_3d, span_3d]/2)

daspect([1 1 1])

xlabel('X [km]')
ylabel('Y [km]')
zlabel('Z [km]')
title(sprintf('Trajectory Approximation (%s Orbit)', Orbit))

legend([h1 h2], {'ODE45', 'LSA'}, 'Location', 'best')

filename = fullfile(results_folder, sprintf('%s_3D_Trajectory_Comparison.png', Orbit));
exportgraphics(gcf, filename, 'Resolution', 300);

% ============================================================
% Position error over normalized orbital period
% ============================================================

figure;
plot(tau_plot, pos_error_km, 'b', 'LineWidth', 1.5)
grid on
set(gcf,'color','w')

xlabel('Normalized Orbit Period')
ylabel('Position Error [km]')
title(sprintf('Position Error: LSA vs ODE45 (%s Orbit)', Orbit))

xlim([0 1])

filename = fullfile(results_folder, sprintf('%s_Position_Error_vs_Period.png', Orbit));
exportgraphics(gcf, filename, 'Resolution', 300);

% ============================================================
% Helper Function
% ============================================================

function DrawMoonCR3BP(MU, erad, prad)

    lstar = 384400;
    npanels = 180;
    alpha = 1;

    folder = 'imagePlanets';
    image_file = imread(fullfile(folder,'moonrev.png'));

    [x, y, z] = ellipsoid((1-MU)*lstar, 0, 0, erad, erad, prad, npanels);
    globe1 = surf(x, y, -z, 'FaceColor', 'none', 'EdgeColor', 0.5*[1 1 1]);

    set(globe1, 'FaceColor', 'texturemap', 'CData', image_file, ...
        'FaceAlpha', alpha, 'EdgeColor', 'none');
end