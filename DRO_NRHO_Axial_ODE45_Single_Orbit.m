clc; 
clear; 
close all;

%********************* Controls ***************************
tspan_type = "open";    % Options "open" or "linspace"
                        % "open" - lets ODE pick the timesteps, may not allow perfect equal spaced in time arcs
                        % "linspace" - uses equal spaced timesteps, should perfectly break into equal space in time arcs

Orbit = "NRHO";          %Options "DRO", "NRHO", or "Axial"

N_revs = 15;

%********************* Parameters ***************************
MU = 0.0121505856;                  % Mass parameter of the Earth-Moon system
lstar = 384400;                     % [km] characteristic length
tstar = 3.751892968837575e+05;      % [s] characteristic time


%********************* Initial Conditions ***************************

% Select timespan type
if strcmp('DRO', Orbit) == 1
    s0 = [1.17;  0;  0;  0; -0.489780292125578;  0];
    time = 3.042534506901380; % best. earlier brian had 3.05
    Color = [0.3010 0.7450 0.9330]; %cyan
    period_colors = flipud(winter(max(N_revs,1)));
    %period_colors = cool(max(N_revs,1));

elseif strcmp('NRHO', Orbit) == 1
    s0= [1.0230033111727; 7.92508194142745e-22; -0.182765473770332; 6.587588408509e-7; -0.10538556046211; -2.449710745305e-6];
    time = 1.524056211242248; %best. earlier brian had 1.53;
    Color = [0.8500 0.3250 0.0980]; % Orange
    period_colors = flipud(autumn(max(N_revs,1)));

elseif strcmp('Axial', Orbit) == 1
    s0 = [0.563868194241744; 0.480452325877622; 5.31556664650877e-30; -0.0914655435406435; -0.100200042489742; -1.05266837879105];
    time = 6.284367749984344; %best
    Color = [0.4660 0.6740 0.1880]; %Green
    period_colors = flipud(parula(max(N_revs,1)));

else
    fprintf("Incorrect Orbit Type")
end


%********************* Reference Trajectory (ODE45) *********************

total_time = N_revs * time;

% Select timespan type
if strcmp('open', tspan_type) == 1
    tspan = [0 total_time];             %May not perfectly split arcs in equal time

elseif strcmp('linspace', tspan_type) == 1
    tspan = linspace(0, total_time, 3000 * N_revs);     %Will perfectly split arcs in equal time

end

% ODE45
options = odeset('RelTol', 1e-13, 'AbsTol', 1e-13); 
[t, s] = ode45(@(t,x) CR3BP(t, x, MU, 1), tspan, s0, options);

% optional: add acceleration (for LCA). only consider [s] for RLCA
for ii = 1:length(t)
    si = s(ii,:);
    s_acc(ii,:) = AccCR3BP(si, MU, 1);
end

Reference = [t, s, s_acc]; %[t, x, y, z, vx, vy, vz, ax, ay, az]


%*********************** Plotting *************************

x_km = s(:,1) * lstar;
y_km = s(:,2) * lstar;
z_km = s(:,3) * lstar;

% Centers
cx = (min(x_km) + max(x_km)) / 2;
cy = (min(y_km) + max(y_km)) / 2;
cz = (min(z_km) + max(z_km)) / 2;

% Spans
span_x = max(x_km) - min(x_km);
span_y = max(y_km) - min(y_km);
span_z = max(z_km) - min(z_km);

pad = 1.10;
span_2d = max(span_x, span_y) * pad;
span_3d = max([span_x, span_y, span_z]) * pad;

min_span = 1e4;
span_2d = max(span_2d, min_span);
span_3d = max(span_3d, min_span);

% Dynamic title
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
    plot3(x_km, y_km, z_km, ...
        'Color', Color, 'LineWidth', 1.5)
else
    for kk = 1:N_revs
        idx = t >= (kk-1)*time & t <= kk*time;
        plot3(x_km(idx), y_km(idx), z_km(idx), ...
            'Color', period_colors(kk,:), 'LineWidth', 1.0)
    end
end

plot3(x_km(1), y_km(1), z_km(1), ...
    'ko', 'MarkerSize', 8, 'LineWidth', 1.5)

plot3(x_km(end), y_km(end), z_km(end), ...
    'rx', 'MarkerSize', 10, 'LineWidth', 1.5)

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
title(plot_title)

if N_revs == 1
    legend('Moon','Trajectory','Initial state','Final state')
else
    colormap(turbo)
    legend([hInit hFinal], {'Initial state','Final state'})
end


%*********************** 3D Plot *************************
figure(2)
DrawMoonCR3BP(MU, 1741, 1741)
hold on

if N_revs == 1
    plot3(x_km, y_km, z_km, ...
        'Color', Color, 'LineWidth', 1.5)
else
    for kk = 1:N_revs
        idx = t >= (kk-1)*time & t <= kk*time;
        plot3(x_km(idx), y_km(idx), z_km(idx), ...
            'Color', period_colors(kk,:), 'LineWidth', 1.0)
    end
end

plot3(x_km(1), y_km(1), z_km(1), ...
    'ko', 'MarkerSize', 8, 'LineWidth', 1.5)

plot3(x_km(end), y_km(end), z_km(end), ...
    'rx', 'MarkerSize', 10, 'LineWidth', 1.5)

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
title(plot_title)

if N_revs == 1
    legend('Moon','Trajectory','Initial state','Final state')
else
    colormap(turbo)
    legend([hInit hFinal], {'Initial state','Final state'})
end


%*********************** Jacobi Constant *************************
JC = CR3BP_JC(s(end,:), MU);
disp(['Jacobi Constant: ', num2str(JC)])


%*********************** Closure Check *************************
initial_state = s(1,:);
final_state   = s(end,:);

state_error_nd = final_state - initial_state;
pos_error_km = norm(state_error_nd(1:3)) * lstar;
vel_error_kms = norm(state_error_nd(4:6)) * (lstar/tstar);

fprintf('\nClosure check for %s:\n', Orbit);
fprintf('Position closure error: %.6e km\n', pos_error_km);
fprintf('Velocity closure error: %.6e km/s\n', vel_error_kms);
fprintf('State error norm (nondimensional): %.6e\n', norm(state_error_nd));


%*********************** Save Dataset *************************

ODE45_Dataset = Reference;

if N_revs == 1
    filename = sprintf('datasets/Single_Orbit_%s_Dataset.mat', Orbit);
else
    filename = sprintf('datasets/Multi_Period_Orbit_%s_Dataset.mat', Orbit);
end

save(filename, 'ODE45_Dataset')

disp(['Workspace cleared, ', filename, ' saved.'])
clearvars -except ODE45_Dataset


% =================================================
% Helper Functions
% =================================================

function ds = CR3BP(t, s, MU, n)
    % Input:
    %   t - Epoch of the integration
    %   s - State (position and velocity only) at the previous epoch
    %   MU - Mass parameter of CR3BP system
    %   n - Mean-motion of the smaller primary around the larger primary,
    %        always 1.0 in any CR3BP
    % Output:
    %   ds - State derivatives (velocity and acc) at the input epoch t, 
    %           using equations-of-motion for CR3BP
    
    % Initialize state at input epoch t
    ds = zeros(6,1);
    
    % Retrieve components from current state
    x = s(1);
    y = s(2);
    z = s(3);
    vx = s(4);
    vy = s(5);
    vz = s(6);
    
    % Compute the distance from the smaller primary
    r = sqrt((x-1+MU)^2+y^2+z^2);

    % Compute the distance from the larger primary
    d = sqrt((x+MU)^2+y^2+z^2); 

    % Assign velocities along x, y and z axes of the CR3BP rotating-frame
    ds(1) = vx; %x-dot
    ds(2) = vy; %y-dot
    ds(3) = vz; %z-dot
    
    % Assign accelerations along x, y and z axes of the CR3BP rotating-frame,
    % using equations of motion for CR3BP model
    ds(4) = (2*n*vy) + x - ((MU*(x-1+MU))/(r^3)) - (((1-MU)*(x+MU))/(d^3)); %x-dot-dot
    ds(5) = (-2*n*vx) + y - ((y*MU)/(r^3)) - ((y*(1-MU))/(d^3)); %y-dot-dot
    ds(6) = -((z*MU)/(r^3))-((z*(1-MU))/(d^3)); %z-dot-dot
    
    %In total ds = [x-dot y-dot z-dot x-dot-dot y-dot-dot z-dot-dot]
end

function Acc = AccCR3BP(s, MU, n)
    % Input:
    %   t - Epoch of the integration
    %   s - State (position and velocity only) at the previous epoch
    %   MU - Mass parameter of CR3BP system
    %   n - Mean-motion of the smaller primary around the larger primary,
    %        always 1.0 in any CR3BP
    % Output:
    %   ds - State derivatives (velocity and acc) at the input epoch t, 
    %           using equations-of-motion for CR3BP
    
    % Initialize state at input epoch t
 
    
    % Retrieve components from current state
    x = s(1);
    y = s(2);
    z = s(3);
    vx = s(4);
    vy = s(5);
    vz = s(6);
    
    % Compute the distance from the smaller primary
    r = sqrt((x-1+MU)^2+y^2+z^2);

    % Compute the distance from the larger primary
    d = sqrt((x+MU)^2+y^2+z^2); 
    
    % Assign accelerations along x, y and z axes of the CR3BP rotating-frame,
    % using equations of motion for CR3BP model
    Acc(1,1) = (2*n*vy) + x - ((MU*(x-1+MU))/(r^3)) - (((1-MU)*(x+MU))/(d^3)); %x-dot-dot
    Acc(1,2) = (-2*n*vx) + y - ((y*MU)/(r^3)) - ((y*(1-MU))/(d^3)); %y-dot-dot
    Acc(1,3) = -((z*MU)/(r^3))-((z*(1-MU))/(d^3)); %z-dot-dot
    
    %In total ds = [x-dot y-dot z-dot x-dot-dot y-dot-dot z-dot-dot]
end

function JC = CR3BP_JC(State, MU)
    x = State(1);
    y = State(2);
    z = State(3);
    xdot = State(4);
    ydot = State(5);
    zdot = State(6);


    % Compute the distance from the smaller primary
    r = sqrt((x-1+MU)^2+y^2+z^2);
    
    % Compute the distance from the larger primary
    d = sqrt((x+MU)^2+y^2+z^2);

    % Psuedo-potential
    U = (1-MU)/d + MU/r + .5*(x^2 + y^2);

    % Jacobi-Constant
    JC = 2*U -(xdot^2 + ydot^2 +zdot^2);
end

% =================================================
% Utilities to draw moon, earth, and trajectories
% =================================================

function DrawMoonCR3BP(MU, erad, prad)
lstar = 384400; % [km] characteristic length

%This code is fulling the image file from a folder in the same directory
%named "imagePlanets"

npanels = 180;   % Number of globe panels around the equator deg/panel = 360/npanels
alpha   = 1; % globe transparency level, 1 = opaque, through 0 = invisible

%image_file = 'https://www.solarsystemscope.com/textures/download/2k_moon.jpg';
folder='imagePlanets';
image_file=imread(fullfile(folder,'moonrev.png'));

[x, y, z] = ellipsoid((1-MU)*lstar, 0, 0, erad, erad, prad, npanels);
globe1 = surf(x, y, -z, 'FaceColor', 'none', 'EdgeColor', 0.5*[1 1 1]);
%cdata = imread(image_file);
set(globe1, 'FaceColor', 'texturemap', 'CData', image_file, 'FaceAlpha', alpha, 'EdgeColor', 'none');

end
