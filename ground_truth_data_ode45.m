clc;
clear;
close all;

% ============================================================
% Periodically Reinitialized Multi-Period Propagation
% Each period starts from the same initial condition s0.
% Global time increases, but the trajectory state resets every period.

% [time, period_id, x, y, z, vx, vy, vz, ax, ay, az]

% Ground truth dataset is generated based on ODE45

% ============================================================

%********************* Controls ***************************
tspan_type = "open";     % Options: "open" or "linspace"
Orbit = "NRHO";           % Options: "DRO", "NRHO", or "Axial"
N_revs = 15;            % Always greater than 1 for this file

%********************* Parameters ***************************
MU = 0.0121505856;
lstar = 384400;                   
tstar = 3.751892968837575e+05;    

%********************* Initial Conditions ***************************
if strcmp('DRO', Orbit)
    s0 = [1.17; 0; 0; 0; -0.489780292125578; 0];
    time = 3.042534506901380;
    Color = [0.3010 0.7450 0.9330];
    period_colors = flipud(winter(N_revs));

elseif strcmp('NRHO', Orbit)
    s0 = [1.0230033111727; 7.92508194142745e-22; -0.182765473770332; ...
          6.587588408509e-7; -0.10538556046211; -2.449710745305e-6];
    time = 1.524056211242248;
    Color = [0.8500 0.3250 0.0980];
    period_colors = flipud(autumn(N_revs));

elseif strcmp('Axial', Orbit)
    s0 = [0.563868194241744; 0.480452325877622; 5.31556664650877e-30; ...
         -0.0914655435406435; -0.100200042489742; -1.05266837879105];
    time = 6.284367749984344;
    Color = [0.4660 0.6740 0.1880];
    period_colors = flipud(parula(N_revs));

else
    error("Incorrect Orbit Type");
end

%********************* Reset Multi-Period Reference Trajectory *********************
options = odeset('RelTol', 1e-13, 'AbsTol', 1e-13);

t_all = [];
s_all = [];
period_id = [];

for kk = 1:N_revs

    if strcmp('open', tspan_type)
        tspan_local = [0 time];

    elseif strcmp('linspace', tspan_type)
        tspan_local = linspace(0, time, 3000);

    else
        error("Incorrect tspan_type");
    end

    % Reinitialize at the same state every period
    [t_local, s_local] = ode45(@(t,x) CR3BP(t, x, MU, 1), ...
                               tspan_local, s0, options);

    % Shift time forward while keeping the state reset
    t_shifted = t_local + (kk - 1)*time;

    % Avoid duplicate time at reset boundaries
    if kk > 1
        t_shifted = t_shifted(2:end);
        s_local = s_local(2:end,:);
    end

    t_all = [t_all; t_shifted];
    s_all = [s_all; s_local];
    period_id = [period_id; kk*ones(size(t_shifted))];

end

t = t_all;
s = s_all;

%********************* Acceleration and Dataset ***************************
s_acc = zeros(length(t), 3);

for ii = 1:length(t)
    s_acc(ii,:) = AccCR3BP(s(ii,:), MU, 1);
end

Reference = [t, period_id, s, s_acc];
% columns: [time, period_id, x, y, z, vx, vy, vz, ax, ay, az]

%*********************** Plotting Preprocessing *************************
x_km = s(:,1) * lstar;
y_km = s(:,2) * lstar;
z_km = s(:,3) * lstar;

cx = (min(x_km) + max(x_km)) / 2;
cy = (min(y_km) + max(y_km)) / 2;
cz = (min(z_km) + max(z_km)) / 2;

span_x = max(x_km) - min(x_km);
span_y = max(y_km) - min(y_km);
span_z = max(z_km) - min(z_km);

pad = 1.10;
span_2d = max(span_x, span_y) * pad;
span_3d = max([span_x, span_y, span_z]) * pad;

min_span = 1e4;
span_2d = max(span_2d, min_span);
span_3d = max(span_3d, min_span);

plot_title = sprintf('%s Repeated Single-Period Propagation (%d Resets)', Orbit, N_revs);

%*********************** 2D Plot *************************
figure(1)
DrawMoonCR3BP(MU, 1741, 1741)
hold on

for kk = 1:N_revs
    idx = period_id == kk;
    plot3(x_km(idx), y_km(idx), z_km(idx), ...
        'Color', period_colors(kk,:), 'LineWidth', 1.0)
end

hInit = plot3(x_km(1), y_km(1), z_km(1), ...
    'ko', 'MarkerSize', 8, 'LineWidth', 1.5);

hFinal = plot3(x_km(end), y_km(end), z_km(end), ...
    'rx', 'MarkerSize', 10, 'LineWidth', 1.5);

set(gcf,'color','w')
axis equal
grid on
view(2)

xlim(cx + [-span_2d, span_2d]/2)
ylim(cy + [-span_2d, span_2d]/2)

xlabel('X [km]')
ylabel('Y [km]')
zlabel('Z [km]')
title([plot_title, ' - 2D View'])
legend([hInit hFinal], {'Initial state','Final state'})

%*********************** 3D Plot *************************
figure(2)
DrawMoonCR3BP(MU, 1741, 1741)
hold on

for kk = 1:N_revs
    idx = period_id == kk;
    plot3(x_km(idx), y_km(idx), z_km(idx), ...
        'Color', period_colors(kk,:), 'LineWidth', 1.0)
end

hInit = plot3(x_km(1), y_km(1), z_km(1), ...
    'ko', 'MarkerSize', 8, 'LineWidth', 1.5);

hFinal = plot3(x_km(end), y_km(end), z_km(end), ...
    'rx', 'MarkerSize', 10, 'LineWidth', 1.5);

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
title([plot_title, ' - 3D View'])
legend([hInit hFinal], {'Initial state','Final state'})

%*********************** Jacobi Constant *************************
JC = CR3BP_JC(s0, MU);
disp(['Jacobi Constant: ', num2str(JC)])

%*********************** Single-Period Closure Check *************************
[t_check, s_check] = ode45(@(t,x) CR3BP(t, x, MU, 1), [0 time], s0, options);

initial_state = s_check(1,:);
final_state   = s_check(end,:);

state_error_nd = final_state - initial_state;
pos_error_km = norm(state_error_nd(1:3)) * lstar;
vel_error_kms = norm(state_error_nd(4:6)) * (lstar/tstar);

fprintf('\nSingle-period closure check for %s:\n', Orbit);
fprintf('Position closure error: %.6e km\n', pos_error_km);
fprintf('Velocity closure error: %.6e km/s\n', vel_error_kms);
fprintf('State error norm (nondimensional): %.6e\n', norm(state_error_nd));

%*********************** Save Dataset *************************
ODE45_Dataset = Reference;

if ~exist('datasets', 'dir')
    mkdir('datasets');
end

filename = sprintf('datasets/GT_%s_%dPeriods_Dataset.mat', Orbit, N_revs);
save(filename, 'ODE45_Dataset')

disp([filename, ' saved.'])

% =================================================
% Helper Functions
% =================================================

function ds = CR3BP(~, s, MU, n)

    ds = zeros(6,1);

    x  = s(1);
    y  = s(2);
    z  = s(3);
    vx = s(4);
    vy = s(5);
    vz = s(6);

    r = sqrt((x - 1 + MU)^2 + y^2 + z^2);
    d = sqrt((x + MU)^2 + y^2 + z^2);

    ds(1) = vx;
    ds(2) = vy;
    ds(3) = vz;

    ds(4) =  2*n*vy + x - (MU*(x - 1 + MU))/(r^3) - ((1 - MU)*(x + MU))/(d^3);
    ds(5) = -2*n*vx + y - (MU*y)/(r^3) - ((1 - MU)*y)/(d^3);
    ds(6) = -(MU*z)/(r^3) - ((1 - MU)*z)/(d^3);
end

function Acc = AccCR3BP(s, MU, n)

    x  = s(1);
    y  = s(2);
    z  = s(3);
    vx = s(4);
    vy = s(5);

    r = sqrt((x - 1 + MU)^2 + y^2 + z^2);
    d = sqrt((x + MU)^2 + y^2 + z^2);

    Acc(1,1) =  2*n*vy + x - (MU*(x - 1 + MU))/(r^3) - ((1 - MU)*(x + MU))/(d^3);
    Acc(1,2) = -2*n*vx + y - (MU*y)/(r^3) - ((1 - MU)*y)/(d^3);
    Acc(1,3) = -(MU*z)/(r^3) - ((1 - MU)*z)/(d^3);
end

function JC = CR3BP_JC(State, MU)

    x = State(1);
    y = State(2);
    z = State(3);
    xdot = State(4);
    ydot = State(5);
    zdot = State(6);

    r = sqrt((x - 1 + MU)^2 + y^2 + z^2);
    d = sqrt((x + MU)^2 + y^2 + z^2);

    U = (1 - MU)/d + MU/r + 0.5*(x^2 + y^2);

    JC = 2*U - (xdot^2 + ydot^2 + zdot^2);
end

function DrawMoonCR3BP(MU, erad, prad)

    lstar = 384400;
    npanels = 180;
    alpha = 1;

    folder = 'imagePlanets';
    image_file = imread(fullfile(folder,'moonrev.png'));

    [x, y, z] = ellipsoid((1 - MU)*lstar, 0, 0, erad, erad, prad, npanels);
    globe1 = surf(x, y, -z, 'FaceColor', 'none', 'EdgeColor', 0.5*[1 1 1]);
    set(globe1, 'FaceColor', 'texturemap', 'CData', image_file, ...
        'FaceAlpha', alpha, 'EdgeColor', 'none');
end