% ============================================================
% single reusable file to generate single orbit dataset for all trajectory
% types: DRO, Axial, NRHO
% can also generate long-term multi-period propagation by setting N_revs
% [t, x, y, z, vx, vy, vz, ax, ay, az]
% ============================================================

clc; 
clear; 
close all;

%********************* Controls ***************************
tspan_type = "open";        % Options: "open" or "linspace"
Solver = "ode45";           % Options: "ode45", "ode78", "ode113"
Orbit = "NRHO";             % Options: "DRO", "NRHO", "Axial"
N_revs = 15;                 % number of periods / revolutions want to be tested

%********************* Parameters ***************************
MU = 0.0121505856;
lstar = 384400;
tstar = 3.751892968837575e+05;

%********************* Initial Conditions ***************************
if strcmp('DRO', Orbit)
    s0 = [1.17; 0; 0; 0; -0.489780292125578; 0];
    time = 3.042534506901380;
    Color = [0.3010 0.7450 0.9330];
    period_colors = flipud(winter(max(N_revs,1)));

elseif strcmp('NRHO', Orbit)
    s0 = [1.0230033111727; 7.92508194142745e-22; -0.182765473770332; ...
          6.587588408509e-7; -0.10538556046211; -2.449710745305e-6];
    time = 1.524056211242248;
    Color = [0.8500 0.3250 0.0980];
    period_colors = flipud(autumn(max(N_revs,1)));

elseif strcmp('Axial', Orbit)
    s0 = [0.563868194241744; 0.480452325877622; 5.31556664650877e-30; ...
         -0.0914655435406435; -0.100200042489742; -1.05266837879105];
    time = 6.284367749984344;
    Color = [0.4660 0.6740 0.1880];
    period_colors = flipud(parula(max(N_revs,1)));

else
    error("Incorrect Orbit Type");
end

%********************* Solver Selection ***************************
switch Solver
    case "ode45"
        ode_solver = @ode45;

    case "ode78"
        if exist('ode78','file') == 2
            ode_solver = @ode78;
        else
            error("ode78 is not found in your MATLAB path. Add your ode78 function file first.");
        end

    case "ode113"
        ode_solver = @ode113;

    otherwise
        error("Incorrect Solver Type. Choose ode45, ode78, or ode113.");
end

%********************* Reference Trajectory ***************************
total_time = N_revs * time;

if strcmp('open', tspan_type)
    tspan = [0 total_time];

elseif strcmp('linspace', tspan_type)
    tspan = linspace(0, total_time, 3000 * N_revs);

else
    error("Incorrect tspan_type");
end

options = odeset('RelTol', 1e-13, 'AbsTol', 1e-13);

fprintf('\nRunning %s for %s over %d period(s)...\n', Solver, Orbit, N_revs);

tic;
[t, s] = ode_solver(@(t,x) CR3BP(t, x, MU, 1), tspan, s0, options);
runtime_seconds = toc;

fprintf('%s runtime: %.6f seconds\n', Solver, runtime_seconds);

%********************* Add Acceleration ***************************
s_acc = zeros(length(t),3);

for ii = 1:length(t)
    s_acc(ii,:) = AccCR3BP(s(ii,:), MU, 1);
end

Reference = [t, s, s_acc]; 
% [t, x, y, z, vx, vy, vz, ax, ay, az]

%*********************** Plotting *************************

set(groot, 'defaultAxesFontName', 'Times New Roman')
set(groot, 'defaultTextFontName', 'Times New Roman')
set(groot, 'defaultLegendFontName', 'Times New Roman')

axisFontSize = 16;
labelFontSize = 16;
titleFontSize = 16;
legendFontSize = 16;

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

if N_revs == 1
    plot_title = sprintf('%s Single Orbit Propagation', Orbit);
else
    plot_title = sprintf('%s Long-Term Propagation (%d Periods)', Orbit, N_revs);
end

%*********************** 2D Plot *************************
figure(1)
DrawMoonCR3BP(MU, 1741, 1741)
hold on

if N_revs == 1
    plot3(x_km, y_km, z_km, 'Color', Color, 'LineWidth', 1.5)
else
    for kk = 1:N_revs
        idx = t >= (kk-1)*time & t <= kk*time;
        plot3(x_km(idx), y_km(idx), z_km(idx), ...
            'Color', period_colors(kk,:), 'LineWidth', 1.0)
    end
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

xlabel('X [km]', 'FontSize', labelFontSize)
ylabel('Y [km]', 'FontSize', labelFontSize)
zlabel('Z [km]', 'FontSize', labelFontSize)

title(plot_title, 'FontSize', titleFontSize)

legend([hInit hFinal], {'Initial state','Final state'}, ...
       'FontSize', legendFontSize)

set(gca, 'FontSize', axisFontSize)

%*********************** 3D Plot *************************
figure(2)
DrawMoonCR3BP(MU, 1741, 1741)
hold on

if N_revs == 1
    plot3(x_km, y_km, z_km, 'Color', Color, 'LineWidth', 1.5)
else
    for kk = 1:N_revs
        idx = t >= (kk-1)*time & t <= kk*time;
        plot3(x_km(idx), y_km(idx), z_km(idx), ...
            'Color', period_colors(kk,:), 'LineWidth', 1.0)
    end
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

xlabel('X [km]', 'FontSize', labelFontSize)
ylabel('Y [km]', 'FontSize', labelFontSize)
zlabel('Z [km]', 'FontSize', labelFontSize)

title(plot_title, 'FontSize', titleFontSize)

legend([hInit hFinal], {'Initial state','Final state'}, ...
       'FontSize', legendFontSize, ...
       'Location', 'best')

set(gca, 'FontSize', axisFontSize)


%*********************** Save Plots *************************
if ~exist('plots', 'dir')
    mkdir('plots');
end

if N_revs == 1
    plotname_2d = sprintf('plots/%s_single_orbit_2D.png', Orbit);
    plotname_3d = sprintf('plots/%s_single_orbit_3D.png', Orbit);
else
    plotname_2d = sprintf('plots/%s_multi_orbit_2D.png', Orbit);
    plotname_3d = sprintf('plots/%s_multi_orbit_3D.png', Orbit);
end

% Tighten 2D figure
figure(1)
set(gca, 'LooseInset', max(get(gca, 'TightInset'), 0.02))
exportgraphics(gca, plotname_2d, 'Resolution', 300)

% Tighten 3D figure
figure(2)
set(gca, 'LooseInset', max(get(gca, 'TightInset'), 0.02))
exportgraphics(gca, plotname_3d, 'Resolution', 300)

fprintf('\nSaved plots:\n%s\n%s\n', plotname_2d, plotname_3d);


%*********************** Jacobi Constant *************************
JC = CR3BP_JC(s(end,:), MU);
disp(['Jacobi Constant: ', num2str(JC)])

%*********************** Closure Check *************************
initial_state = s(1,:);
final_state   = s(end,:);

state_error_nd = final_state - initial_state;
pos_error_km = norm(state_error_nd(1:3)) * lstar;
vel_error_kms = norm(state_error_nd(4:6)) * (lstar/tstar);

fprintf('\nClosure check for %s using %s:\n', Orbit, Solver);
fprintf('Position closure error: %.6e km\n', pos_error_km);
fprintf('Velocity closure error: %.6e km/s\n', vel_error_kms);
fprintf('State error norm (nondimensional): %.6e\n', norm(state_error_nd));

%*********************** Save Dataset *************************
Dataset = Reference;

if ~exist('datasets', 'dir')
    mkdir('datasets');
end

if N_revs == 1
    filename = sprintf('datasets/%s_Single_Orbit_%s_Dataset.mat', Solver, Orbit);
else
    filename = sprintf('datasets/%s_Multi_Period_%s_%dPeriods_Dataset.mat', Solver, Orbit, N_revs);
end

save(filename, 'Dataset', 'runtime_seconds', 'Solver', 'Orbit', 'N_revs')

fprintf('\nSaved dataset:\n%s\n', filename);

clearvars -except Dataset runtime_seconds Solver Orbit N_revs

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

    r = sqrt((x-1+MU)^2 + y^2 + z^2);
    d = sqrt((x+MU)^2 + y^2 + z^2); 

    ds(1) = vx;
    ds(2) = vy;
    ds(3) = vz;

    ds(4) =  2*n*vy + x - (MU*(x-1+MU))/(r^3) - ((1-MU)*(x+MU))/(d^3);
    ds(5) = -2*n*vx + y - (MU*y)/(r^3) - ((1-MU)*y)/(d^3);
    ds(6) = -(MU*z)/(r^3) - ((1-MU)*z)/(d^3);
end

function Acc = AccCR3BP(s, MU, n)

    x  = s(1);
    y  = s(2);
    z  = s(3);
    vx = s(4);
    vy = s(5);

    r = sqrt((x-1+MU)^2 + y^2 + z^2);
    d = sqrt((x+MU)^2 + y^2 + z^2); 
    
    Acc(1,1) =  2*n*vy + x - (MU*(x-1+MU))/(r^3) - ((1-MU)*(x+MU))/(d^3);
    Acc(1,2) = -2*n*vx + y - (MU*y)/(r^3) - ((1-MU)*y)/(d^3);
    Acc(1,3) = -(MU*z)/(r^3) - ((1-MU)*z)/(d^3);
end

function JC = CR3BP_JC(State, MU)

    x = State(1);
    y = State(2);
    z = State(3);
    xdot = State(4);
    ydot = State(5);
    zdot = State(6);

    r = sqrt((x-1+MU)^2 + y^2 + z^2);
    d = sqrt((x+MU)^2 + y^2 + z^2);

    U = (1-MU)/d + MU/r + 0.5*(x^2 + y^2);

    JC = 2*U - (xdot^2 + ydot^2 + zdot^2);
end

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