clc;
clear;
close all;

% ============================================================
% Multi-period comparison:
% GT vs ODE45 multi-period vs LSA repeated per orbit period
% ============================================================

Orbit = 'Axial';
Dataset_Periods = 100;     % periods in filename
Periods_To_Test = 6;      % periods to actually analyze

gt_file  = sprintf('datasets/GT_%s_%dPeriods_Dataset.mat', Orbit, Dataset_Periods);
ode_file = sprintf('datasets/ode45_Multi_Period_%s_%dPeriods_Dataset.mat', Orbit, Dataset_Periods);
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
% Load GT dataset
% GT format:
% [time, period_id, x, y, z, vx, vy, vz, ax, ay, az]
% ============================================================

GT = load(gt_file);

if isfield(GT, 'ODE45_Dataset')
    GT_Dataset = GT.ODE45_Dataset;
elseif isfield(GT, 'Dataset')
    GT_Dataset = GT.Dataset;
else
    error('GT file does not contain ODE45_Dataset or Dataset.');
end

if isfield(GT, 'tstar')
    tstar = GT.tstar;
end

% Keep only selected number of periods
GT_Dataset = GT_Dataset(GT_Dataset(:,2) <= Periods_To_Test, :);

t_gt = GT_Dataset(:,1);
period_id = GT_Dataset(:,2);
X_gt = GT_Dataset(:,3:11);

pos_gt = X_gt(:,1:3);
vel_gt = X_gt(:,4:6);
acc_gt = X_gt(:,7:9);

% ============================================================
% Load ODE45 multi-period dataset
% Expected format:
% [time, x, y, z, vx, vy, vz, ax, ay, az]
% or
% [time, period_id, x, y, z, vx, vy, vz, ax, ay, az]
% ============================================================

ODE = load(ode_file);

if isfield(ODE, 'Dataset')
    ODE_Dataset = ODE.Dataset;
elseif isfield(ODE, 'ODE45_Dataset')
    ODE_Dataset = ODE.ODE45_Dataset;
else
    error('ODE file does not contain Dataset or ODE45_Dataset.');
end

if size(ODE_Dataset,2) == 10
    t_ode = ODE_Dataset(:,1);
    X_ode_raw = ODE_Dataset(:,2:10);
elseif size(ODE_Dataset,2) == 11
    t_ode = ODE_Dataset(:,1);
    X_ode_raw = ODE_Dataset(:,3:11);
else
    error('Unexpected ODE dataset format.');
end

% Interpolate ODE45 onto GT timestamps
X_ode = interp1(t_ode, X_ode_raw, t_gt, 'linear', 'extrap');

pos_ode = X_ode(:,1:3);
vel_ode = X_ode(:,4:6);
acc_ode = X_ode(:,7:9);

% ============================================================
% Load LSA polynomial coefficients
% ============================================================

LSA = load(lsa_file);
Final_Table = LSA.Final_Table;

components = string(Final_Table.Component);
best_degrees = Final_Table.BestDegree;

expected_components = ["x","y","z","vx","vy","vz","ax","ay","az"];

num_points = length(t_gt);
num_components = 9;

X_lsa = zeros(num_points, num_components);

% ============================================================
% Build orbit-axis and local normalized time per period
% ============================================================

unique_periods = unique(period_id, 'stable');
num_periods_found = length(unique_periods);

orbit_axis = zeros(num_points,1);
tau_lsa = zeros(num_points,1);

for p = 1:num_periods_found

    pid = unique_periods(p);
    idx = period_id == pid;

    t_local = t_gt(idx);

    tau_local = 2*(t_local - t_local(1)) / (t_local(end) - t_local(1)) - 1;
    tau_lsa(idx) = tau_local;

    orbit_frac = (t_local - t_local(1)) / (t_local(end) - t_local(1));
    orbit_axis(idx) = (p-1) + orbit_frac;
end

% ============================================================
% Reconstruct LSA over all GT timestamps
% Polynomial is evaluated locally within each orbit period
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
% LSA vs GT and ODE45 vs GT
% ============================================================

pos_err_lsa_km = sqrt(sum((pos_lsa - pos_gt).^2, 2)) * lstar;
pos_err_ode_km = sqrt(sum((pos_ode - pos_gt).^2, 2)) * lstar;

vel_err_lsa_kms = sqrt(sum((vel_lsa - vel_gt).^2, 2)) * (lstar / tstar);
vel_err_ode_kms = sqrt(sum((vel_ode - vel_gt).^2, 2)) * (lstar / tstar);

acc_err_lsa_kms2 = sqrt(sum((acc_lsa - acc_gt).^2, 2)) * (lstar / tstar^2);
acc_err_ode_kms2 = sqrt(sum((acc_ode - acc_gt).^2, 2)) * (lstar / tstar^2);

fprintf('\n============================================================\n');
fprintf('Multi-period comparison for %s\n', Orbit);
fprintf('Dataset periods : %d\n', Dataset_Periods);
fprintf('Tested periods  : %d\n', num_periods_found);
fprintf('============================================================\n');

fprintf('\nPosition RMSE:\n');
fprintf('LSA  vs GT : %.6e km\n', sqrt(mean(pos_err_lsa_km.^2)));
fprintf('ODE45 vs GT: %.6e km\n', sqrt(mean(pos_err_ode_km.^2)));

fprintf('\nVelocity RMSE:\n');
fprintf('LSA  vs GT : %.6e km/s\n', sqrt(mean(vel_err_lsa_kms.^2)));
fprintf('ODE45 vs GT: %.6e km/s\n', sqrt(mean(vel_err_ode_kms.^2)));

fprintf('\nAcceleration RMSE:\n');
fprintf('LSA  vs GT : %.6e km/s^2\n', sqrt(mean(acc_err_lsa_kms2.^2)));
fprintf('ODE45 vs GT: %.6e km/s^2\n', sqrt(mean(acc_err_ode_kms2.^2)));

% ============================================================
% Convert trajectories to km
% ============================================================

x_gt_km = pos_gt(:,1) * lstar;
y_gt_km = pos_gt(:,2) * lstar;
z_gt_km = pos_gt(:,3) * lstar;

x_ode_km = pos_ode(:,1) * lstar;
y_ode_km = pos_ode(:,2) * lstar;
z_ode_km = pos_ode(:,3) * lstar;

x_lsa_km = pos_lsa(:,1) * lstar;
y_lsa_km = pos_lsa(:,2) * lstar;
z_lsa_km = pos_lsa(:,3) * lstar;

all_x = [x_gt_km; x_ode_km; x_lsa_km];
all_y = [y_gt_km; y_ode_km; y_lsa_km];
all_z = [z_gt_km; z_ode_km; z_lsa_km];

cx = (min(all_x) + max(all_x)) / 2;
cy = (min(all_y) + max(all_y)) / 2;
cz = (min(all_z) + max(all_z)) / 2;

span_x = max(all_x) - min(all_x);
span_y = max(all_y) - min(all_y);
span_z = max(all_z) - min(all_z);

pad = 1.10;
span_3d = max([span_x, span_y, span_z]) * pad;
span_3d = max(span_3d, 1e4);

% ============================================================
% 3D trajectory comparison
% ============================================================

figure;
DrawMoonCR3BP(MU, 1741, 1741)
hold on

h1 = plot3(x_ode_km, y_ode_km, z_ode_km, 'b', 'LineWidth', 1.2);
h2 = plot3(x_lsa_km, y_lsa_km, z_lsa_km, 'r', 'LineWidth', 1.2);
h3 = plot3(x_gt_km, y_gt_km, z_gt_km, 'g', 'LineWidth', 1.2);

plot3(x_gt_km(1), y_gt_km(1), z_gt_km(1), 'ko', 'MarkerSize', 8, 'LineWidth', 1.5);
plot3(x_gt_km(end), y_gt_km(end), z_gt_km(end), 'kx', 'MarkerSize', 10, 'LineWidth', 1.5);

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
title(sprintf('%s Multi-period Trajectory Comparison (%d Periods)', Orbit, num_periods_found))

legend([h1 h2 h3], {'ODE45', 'LSA', 'GT'}, 'Location', 'best')

filename = fullfile(results_folder, sprintf('%s_%dPeriods_3D_Trajectory_Comparison.png', Orbit, num_periods_found));
exportgraphics(gcf, filename, 'Resolution', 300);

% ============================================================
% Position error over orbit periods
% ============================================================

figure;
hold on
grid on

plot(orbit_axis, pos_err_ode_km, 'b', 'LineWidth', 1.4)
plot(orbit_axis, pos_err_lsa_km, 'r', 'LineWidth', 1.4)

for p = 1:num_periods_found
    xline(p, ':k', 'HandleVisibility', 'off');
end

set(gcf,'color','w')

xlabel('Orbit Period')
ylabel('Position Error [km]')
title(sprintf('%s Position Error over %d Periods', Orbit, num_periods_found))

xlim([0 num_periods_found])
xticks(0:num_periods_found)

legend({'ODE45 vs GT', 'LSA vs GT'}, 'Location', 'best')

filename = fullfile(results_folder, sprintf('%s_%dPeriods_Position_Error_vs_Orbit.png', Orbit, num_periods_found));
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