clc;
clear;
close all;

% Load MATLAB file
orbit = 'Axial';
degrees = [6, 10, 13];

mat = load(sprintf('datasets/ode45_Single_Orbit_%s_Dataset.mat', orbit));

data = double(mat.Dataset);

t = data(:, 1);
data = data(:, 2:end);

N = size(data, 1);

MU = 0.0121505856;
lstar = 384400;
tstar = 3.751892968837575e+05;

state_data = data;
col_names = {'x', 'y', 'z', 'vx', 'vy', 'vz', 'ax', 'ay', 'az'};

x_ref_km = data(:,1) * lstar;
y_ref_km = data(:,2) * lstar;
z_ref_km = data(:,3) * lstar;

cx = (min(x_ref_km) + max(x_ref_km)) / 2;
cy = (min(y_ref_km) + max(y_ref_km)) / 2;
cz = (min(z_ref_km) + max(z_ref_km)) / 2;

span_x = max(x_ref_km) - min(x_ref_km);
span_y = max(y_ref_km) - min(y_ref_km);
span_z = max(z_ref_km) - min(z_ref_km);

pad = 1.10;
span_3d = max([span_x, span_y, span_z]) * pad;

min_span = 1e4;
span_3d = max(span_3d, min_span);

figure;
set(gcf,'color','w')
tiledlayout(1,3);

for d = 1:length(degrees)

    degree = degrees(d);

    A = fliplr(vander(t));
    A = A(:, 1:degree+1);

    ATA = A.' * A;
    ATA_inv = inv(ATA);
    ATA_inv_AT = ATA_inv * A.';

    prediction = {};
    fprintf("\nManual Least-Squares Results: Degree %d\n\n", degree);

    for i = 1:size(state_data, 2)
        b = state_data(:, i);

        coefficients = ATA_inv_AT * b;
        pred = A * coefficients;
        prediction{end+1} = pred;

        residual = b - pred;
        rmse = sqrt(mean(residual.^2));

        fprintf('%3s | RMSE = %.2e NDU\n', col_names{i}, rmse);
    end

    prediction = cell2mat(prediction);

    x_fit_km = prediction(:,1) * lstar;
    y_fit_km = prediction(:,2) * lstar;
    z_fit_km = prediction(:,3) * lstar;

    nexttile;
    DrawMoonCR3BP(MU, 1741, 1741)
    hold on

    plot3(x_ref_km, y_ref_km, z_ref_km, ...
        'Color', [0.2745 0.5098 0.7059], ...
        'LineWidth', 1.8, ...
        'DisplayName', 'ODE45 Reference');

    plot3(x_fit_km, y_fit_km, z_fit_km, ...
        'r', ...
        'LineWidth', 1.8, ...
        'DisplayName', 'LSA Fit');

    plot3(x_ref_km(1), y_ref_km(1), z_ref_km(1), ...
        'ko', 'MarkerSize', 8, 'LineWidth', 1.5, ...
        'DisplayName', 'Initial state');

    plot3(x_fit_km(end), y_fit_km(end), z_fit_km(end), ...
        'rx', 'MarkerSize', 10, 'LineWidth', 1.5, ...
        'DisplayName', 'Final fitted state');

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
    title(sprintf('Degree %d', degree))
    legend('Location','southoutside')
end

function DrawMoonCR3BP(MU, erad, prad)

    lstar = 384400;
    npanels = 180;
    alpha = 1;

    folder = 'imagePlanets';
    image_file = imread(fullfile(folder,'moonrev.png'));

    [x, y, z] = ellipsoid((1-MU)*lstar, 0, 0, erad, erad, prad, npanels);

    globe1 = surf(x, y, -z, ...
        'FaceColor', 'none', ...
        'EdgeColor', 0.5*[1 1 1]);

    set(globe1, ...
        'FaceColor', 'texturemap', ...
        'CData', image_file, ...
        'FaceAlpha', alpha, ...
        'EdgeColor', 'none');

    % remove moon from legend
    globe1.Annotation.LegendInformation.IconDisplayStyle = 'off';
end