clc;
clear;
close all;

% ============================================================
% Closure-period scan for CR3BP trajectories
% ============================================================

Orbit = "NRHO";   % Options: "DRO", "NRHO", "Axial"

MU = 0.0121505856;
lstar = 384400;                   
tstar = 3.751892968837575e+05;    

% ---------------- Initial conditions and scan range ----------
switch Orbit
    case "DRO"
        s0 = [1.17; 0; 0; 0; -0.489780292125578; 0];
        t_min = 3.04;
        t_max = 3.045;

    case "NRHO"
        s0 = [1.0230033111727; 7.92508194142745e-22; -0.182765473770332; ...
              6.587588408509e-7; -0.10538556046211; -2.449710745305e-6];
        t_min = 1.523;
        t_max = 1.525;

    case "Axial"
        s0 = [0.563868194241744; 0.480452325877622; 5.31556664650877e-30; ...
             -0.0914655435406435; -0.100200042489742; -1.05266837879105];
        t_min = 6.284;
        t_max = 6.285;

    otherwise
        error("Incorrect Orbit Type");
end

% ---------------- Scan settings -----------------------------
N_scan = 5000;  
time_values = linspace(t_min, t_max, N_scan);

options = odeset('RelTol', 1e-13, 'AbsTol', 1e-13);

pos_error_km = zeros(N_scan,1);
vel_error_kms = zeros(N_scan,1);
state_error_nd = zeros(N_scan,1);

% ---------------- Period scan -------------------------------
for i = 1:N_scan
    tf = time_values(i);

    [~, s] = ode45(@(t,x) CR3BP(t, x, MU, 1), [0 tf], s0, options);

    err = s(end,:).' - s0;

    pos_error_km(i) = norm(err(1:3)) * lstar;
    vel_error_kms(i) = norm(err(4:6)) * (lstar/tstar);
    state_error_nd(i) = norm(err);
end

% ---------------- Best closure result ------------------------
[best_error, best_idx] = min(state_error_nd);
best_time = time_values(best_idx);

fprintf('\nBest closure result for %s:\n', Orbit);
fprintf('Best time: %.15f nondimensional time units\n', best_time);
fprintf('Position closure error: %.6e km\n', pos_error_km(best_idx));
fprintf('Velocity closure error: %.6e km/s\n', vel_error_kms(best_idx));
fprintf('State error norm: %.6e\n', best_error);

% Show top 10 candidates
results = table(time_values.', pos_error_km, vel_error_kms, state_error_nd, ...
    'VariableNames', {'Time', 'PositionError_km', 'VelocityError_km_s', 'StateError_ND'});

results = sortrows(results, 'StateError_ND');

disp(' ');
disp('Top 10 candidate closure times:');
disp(results(1:10,:));

% ---------------- Plot error vs candidate period --------------
figure;
plot(time_values, pos_error_km, 'LineWidth', 1.5);
hold on;
plot(best_time, pos_error_km(best_idx), 'ro', 'MarkerSize', 8, 'LineWidth', 1.5);
grid on;
xlabel('Final propagation time');
ylabel('Position closure error [km]');
title(sprintf('%s closure scan: position error', Orbit));

figure;
plot(time_values, state_error_nd, 'LineWidth', 1.5);
hold on;
plot(best_time, best_error, 'ro', 'MarkerSize', 8, 'LineWidth', 1.5);
grid on;
xlabel('Final propagation time');
ylabel('Full state closure error [nondimensional]');
title(sprintf('%s closure scan: full state error', Orbit));

% ============================================================
% Helper function
% ============================================================

function ds = CR3BP(~, s, MU, n)

    ds = zeros(6,1);

    x  = s(1);
    y  = s(2);
    z  = s(3);
    vx = s(4);
    vy = s(5);

    r = sqrt((x - 1 + MU)^2 + y^2 + z^2);  % distance to Moon
    d = sqrt((x + MU)^2 + y^2 + z^2);      % distance to Earth

    ds(1) = s(4);
    ds(2) = s(5);
    ds(3) = s(6);

    ds(4) =  2*n*vy + x - (MU*(x - 1 + MU))/(r^3) - ((1 - MU)*(x + MU))/(d^3);
    ds(5) = -2*n*vx + y - (MU*y)/(r^3) - ((1 - MU)*y)/(d^3);
    ds(6) = -(MU*z)/(r^3) - ((1 - MU)*z)/(d^3);
end