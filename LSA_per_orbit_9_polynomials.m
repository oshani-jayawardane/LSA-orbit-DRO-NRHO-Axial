% ============================================================
% LSA is fitted using the dataset generated using ODE45
% ============================================================

clc;
clear;
close all;

orbits = ["DRO","NRHO","Axial"];

max_degree = 20;
K = 5;
cond_threshold = 1e12;

component_names = ["x","y","z","vx","vy","vz","ax","ay","az"];

% Dimensional constants
lstar = 384400;                         % km
default_tstar = 3.751892968837575e+05;  % s

if ~exist('LSA', 'dir')
    mkdir('LSA');
end

% Summary storage
RMSE_Summary = table( ...
    strings(length(orbits),1), ...
    zeros(length(orbits),1), zeros(length(orbits),1), zeros(length(orbits),1), ...
    zeros(length(orbits),1), zeros(length(orbits),1), zeros(length(orbits),1), ...
    'VariableNames', { ...
    'Orbit', ...
    'Position_RMSE_NDU', 'Velocity_RMSE_NDU', 'Acceleration_RMSE_NDU', ...
    'Position_RMSE_km', 'Velocity_RMSE_km_s', 'Acceleration_RMSE_km_s2'});

figure;
tiledlayout(3,1);

for o = 1:length(orbits)

    Orbit = orbits(o);

    mat_file = sprintf("datasets/ode45_Single_Orbit_%s_Dataset.mat", Orbit);
    output_file_mat = sprintf("LSA/LSA_%s_Table.mat", Orbit);
    output_file_csv = sprintf("LSA/LSA_%s_Table.csv", Orbit);

    S = load(mat_file);

    % ========================================================
    % Load Dataset
    % Expected format: [t, x, y, z, vx, vy, vz, ax, ay, az]
    % ========================================================

    data = S.Dataset;

    if isfield(S, 'tstar')
        tstar = S.tstar;
    else
        tstar = default_tstar;
    end

    t_raw = data(:,1);
    Y = data(:,2:end);

    t_min = min(t_raw);
    t_max = max(t_raw);
    t = 2*(t_raw - t_min)/(t_max - t_min) - 1;

    n_samples = length(t);
    n_components = size(Y,2);

    best_degrees = zeros(n_components,1);
    coefficients = NaN(n_components,max_degree+1);
    component_rmse_nd = zeros(n_components,1);
    component_rmse_dim = zeros(n_components,1);

    Y_fit = zeros(size(Y));
    cv_rmse_by_component = NaN(n_components,max_degree+1);

    cv = cvpartition(n_samples,'KFold',K);

    for c = 1:n_components

        y = Y(:,c);
        degree_rmse = NaN(max_degree+1,1);
        ill_conditioned = false(max_degree+1,1);

        for deg = 0:max_degree

            fold_rmse = NaN(K,1);
            degree_bad = false;

            for k = 1:K

                train_idx = training(cv,k);
                val_idx = test(cv,k);

                V_train = buildVandermonde(t(train_idx),deg);
                V_val = buildVandermonde(t(val_idx),deg);

                y_train = y(train_idx);
                y_val = y(val_idx);

                M = V_train' * V_train;

                if cond(M) > cond_threshold || rank(M) < deg+1
                    degree_bad = true;
                    fprintf("Ill-conditioned: Orbit=%s | %s | degree=%d | fold=%d\n", ...
                        Orbit, component_names(c), deg, k);
                end

                a_hat = V_train \ y_train;
                y_pred = V_val * a_hat;

                fold_rmse(k) = sqrt(mean((y_val - y_pred).^2));
            end

            degree_rmse(deg+1) = mean(fold_rmse,'omitnan');
            ill_conditioned(deg+1) = degree_bad;
        end

        valid = ~ill_conditioned;

        if any(valid)
            temp = degree_rmse;
            temp(~valid) = Inf;
            [~, best_idx] = min(temp);
            best_deg = best_idx - 1;
        else
            error("No stable degree found for Orbit=%s, Component=%s", Orbit, component_names(c));
        end

        best_degrees(c) = best_deg;
        cv_rmse_by_component(c,:) = degree_rmse.';

        % ====================================================
        % Fit full dataset using selected best degree
        % ====================================================

        V_final = buildVandermonde(t,best_deg);
        a_final = V_final \ y;
        y_fit = V_final * a_final;

        Y_fit(:,c) = y_fit;
        coefficients(c,1:length(a_final)) = a_final';

        % Component RMSE in nondimensional units
        component_rmse_nd(c) = sqrt(mean((y - y_fit).^2));

        % Component RMSE in dimensional units
        if c <= 3
            component_rmse_dim(c) = component_rmse_nd(c) * lstar;
        elseif c <= 6
            component_rmse_dim(c) = component_rmse_nd(c) * (lstar / tstar);
        else
            component_rmse_dim(c) = component_rmse_nd(c) * (lstar / tstar^2);
        end
    end

    % ========================================================
    % Vector RMSE for position, velocity, acceleration
    % ========================================================

    pos_error = Y(:,1:3) - Y_fit(:,1:3);
    vel_error = Y(:,4:6) - Y_fit(:,4:6);
    acc_error = Y(:,7:9) - Y_fit(:,7:9);

    position_rmse_nd = sqrt(mean(sum(pos_error.^2,2)));
    velocity_rmse_nd = sqrt(mean(sum(vel_error.^2,2)));
    acceleration_rmse_nd = sqrt(mean(sum(acc_error.^2,2)));

    position_rmse_km = position_rmse_nd * lstar;
    velocity_rmse_kms = velocity_rmse_nd * (lstar / tstar);
    acceleration_rmse_kms2 = acceleration_rmse_nd * (lstar / tstar^2);

    % Store summary values
    RMSE_Summary.Orbit(o) = Orbit;
    RMSE_Summary.Position_RMSE_NDU(o) = position_rmse_nd;
    RMSE_Summary.Velocity_RMSE_NDU(o) = velocity_rmse_nd;
    RMSE_Summary.Acceleration_RMSE_NDU(o) = acceleration_rmse_nd;
    RMSE_Summary.Position_RMSE_km(o) = position_rmse_km;
    RMSE_Summary.Velocity_RMSE_km_s(o) = velocity_rmse_kms;
    RMSE_Summary.Acceleration_RMSE_km_s2(o) = acceleration_rmse_kms2;

    fprintf('\n=====================================================\n');
    fprintf('LSA RESULTS FOR %s\n', Orbit);
    fprintf('=====================================================\n');

    fprintf('\nComponent-wise RMSE:\n');
    for c = 1:n_components

        if c <= 3
            unit_str = 'km';
        elseif c <= 6
            unit_str = 'km/s';
        else
            unit_str = 'km/s^2';
        end

        fprintf('%3s | Best Degree = %2d | RMSE NDU = %.6e | RMSE Dim = %.6e %s\n', ...
            component_names(c), best_degrees(c), component_rmse_nd(c), ...
            component_rmse_dim(c), unit_str);
    end

    fprintf('-----------------------------------------------------\n');
    fprintf('Vector RMSE in NDU:\n');
    fprintf('Position RMSE     = %.6e\n', position_rmse_nd);
    fprintf('Velocity RMSE     = %.6e\n', velocity_rmse_nd);
    fprintf('Acceleration RMSE = %.6e\n', acceleration_rmse_nd);

    fprintf('\nVector RMSE in dimensional units:\n');
    fprintf('Position RMSE     = %.6e km\n', position_rmse_km);
    fprintf('Velocity RMSE     = %.6e km/s\n', velocity_rmse_kms);
    fprintf('Acceleration RMSE = %.6e km/s^2\n', acceleration_rmse_kms2);
    fprintf('=====================================================\n');

    % ========================================================
    % Save coefficient table
    % ========================================================

    header = ["Component","BestDegree"];
    for j = 0:max_degree
        header = [header, "c"+string(j)];
    end

    Final_Cell = cell(n_components,length(header));

    for i = 1:n_components
        Final_Cell{i,1} = component_names(i);
        Final_Cell{i,2} = best_degrees(i);

        for j = 1:max_degree+1
            Final_Cell{i,j+2} = coefficients(i,j);
        end
    end

    Final_Table = cell2table(Final_Cell, ...
        'VariableNames', matlab.lang.makeValidName(header));

    save(output_file_mat,"Final_Table");
    writetable(Final_Table,output_file_csv);

    fprintf('\nSaved:\n%s\n%s\n', output_file_mat, output_file_csv);


    % ============================================================
    % Plot settings
    % ============================================================

    set(groot, 'defaultAxesFontName', 'Times New Roman')
    set(groot, 'defaultTextFontName', 'Times New Roman')
    set(groot, 'defaultLegendFontName', 'Times New Roman')
    
    fontName = 'Times New Roman';
    
    axisFontSize = 10;
    labelFontSize = 10;
    titleFontSize = 12;
    legendFontSize = 12;
    sgtitleFontSize = 14;
    
    lineWidth = 1.2;
    markerSize = 6;
    
    set(groot, 'defaultAxesFontName', fontName)
    set(groot, 'defaultTextFontName', fontName)
    set(groot, 'defaultLegendFontName', fontName)

    % ========================================================
    % CV RMSE plot
    % ========================================================

    nexttile;
    hold on;
    grid on;

    degrees = 0:max_degree;

    for c = 1:n_components
        plot(degrees, cv_rmse_by_component(c,:), ...
            'LineWidth', lineWidth, ...
            'DisplayName', component_names(c));
    end
    
    for c = 1:n_components
        bd = best_degrees(c);
        plot(bd, cv_rmse_by_component(c,bd+1), 'o', ...
            'MarkerSize', markerSize, ...
            'LineWidth', lineWidth, ...
            'HandleVisibility','off');
    end
    
    xlabel('Polynomial Degree', 'FontSize', labelFontSize);
    ylabel('Val RMSE [NDU]', 'FontSize', labelFontSize);
    title(sprintf('%s Orbit', Orbit), 'FontSize', titleFontSize);
    
    set(gca, 'FontSize', axisFontSize, 'FontName', fontName);
    
    if o == 1
        legend('Location','eastoutside', ...
               'FontSize', legendFontSize, ...
               'FontName', fontName);
    end
end

% ============================================================
% Save summary RMSE table
% ============================================================

save("LSA/LSA_Vector_RMSE_Summary.mat", "RMSE_Summary");
writetable(RMSE_Summary, "LSA/LSA_Vector_RMSE_Summary.csv");

disp(" ");
disp("Saved RMSE summary table:");
disp(RMSE_Summary);

set(gcf,'color','w');

sgtitle('Best-Fit Polynomial Degree Selection using 5-Fold Cross-Validation', ...
        'FontName', fontName, ...
        'FontSize', sgtitleFontSize);

saveas(gcf, "LSA/All_Orbits_CV_RMSE_vs_Degree.png");

clearvars -except Final_Table RMSE_Summary

function V = buildVandermonde(t,degree)

n = length(t);
V = zeros(n,degree+1);

for k = 0:degree
    V(:,k+1) = t.^k;
end

end