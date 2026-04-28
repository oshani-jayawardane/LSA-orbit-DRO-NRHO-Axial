clc;
clear;
close all;

% ============================================================
% Position RMSE and Condition Number vs Polynomial Degree
% Dataset format: [t, x, y, z, vx, vy, vz, ax, ay, az]
% ============================================================

orbits = {'NRHO', 'DRO', 'Axial'};
degrees = 3:15;

lstar = 384400;

position_rmse_all = zeros(length(orbits), length(degrees));
condition_number_all = zeros(length(orbits), length(degrees));

for o = 1:length(orbits)

    orbit = orbits{o};

    mat = load(sprintf('datasets/ode45_Single_Orbit_%s_Dataset.mat', orbit));

    data = double(mat.Dataset);

    % Use actual time column
    t = data(:,1);

    % Drop time column
    data = data(:,2:end);

    state_data = data;

    fprintf('\nPosition RMSE vs Degree for %s\n\n', orbit);

    for d = 1:length(degrees)

        degree = degrees(d);

        % Vandermonde matrix
        A = fliplr(vander(t));
        A = A(:,1:degree+1);

        % Normal equation terms
        ATA = A.' * A;

        % Condition number of normal matrix
        condition_number_all(o,d) = cond(ATA);

        ATA_inv = inv(ATA);
        ATA_inv_AT = ATA_inv * A.';

        prediction = zeros(size(state_data));

        % Fit all 9 components using same degree
        for i = 1:size(state_data,2)
            b = state_data(:,i);
            coefficients = ATA_inv_AT * b;
            prediction(:,i) = A * coefficients;
        end

        % Position RMSE using x, y, z only
        pos_error = state_data(:,1:3) - prediction(:,1:3);

        position_rmse = sqrt(mean(sum(pos_error.^2,2)));

        position_rmse_all(o,d) = position_rmse;

        fprintf('Degree %2d | Position RMSE = %.6e NDU | cond(ATA) = %.6e\n', ...
            degree, position_rmse, condition_number_all(o,d));
    end
end

% ============================================================
% Create tables
% ============================================================

header = strings(1, length(degrees));
for i = 1:length(degrees)
    header(i) = "Degree_" + string(degrees(i));
end

RMSE_Table = array2table(position_rmse_all, ...
    'VariableNames', matlab.lang.makeValidName(header), ...
    'RowNames', orbits);

Condition_Table = array2table(condition_number_all, ...
    'VariableNames', matlab.lang.makeValidName(header), ...
    'RowNames', orbits);

disp('Position RMSE Table')
disp(RMSE_Table)

disp('Condition Number Table')
disp(Condition_Table)

% ============================================================
% Optional save
% ============================================================

if ~exist('LSA', 'dir')
    mkdir('LSA');
end

rmse_output_file = 'LSA/All_Orbits_Position_RMSE_Degrees_3_to_15.csv';
cond_output_file = 'LSA/All_Orbits_Condition_Number_Degrees_3_to_15.csv';

writetable(RMSE_Table, rmse_output_file, 'WriteRowNames', true);
writetable(Condition_Table, cond_output_file, 'WriteRowNames', true);

fprintf('\nSaved RMSE table to:\n%s\n', rmse_output_file);
fprintf('\nSaved condition number table to:\n%s\n', cond_output_file);