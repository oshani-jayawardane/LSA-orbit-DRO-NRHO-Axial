clc;
clear;
close all;

orbits = ["DRO","NRHO","Axial"];

max_degree = 20;
K = 5;
cond_threshold = 1e12;
data_var = "ODE45_Dataset";
component_names = ["x","y","z","vx","vy","vz","ax","ay","az"];

if ~exist('LSA', 'dir')
    mkdir('LSA');
end

figure;
tiledlayout(3,1);

for o = 1:length(orbits)

    Orbit = orbits(o);

    mat_file = sprintf("datasets/Single_Orbit_%s_Dataset.mat", Orbit);
    output_file_mat = sprintf("LSA/LSA_%s_Table.mat", Orbit);
    output_file_csv = sprintf("LSA/LSA_%s_Table.csv", Orbit);

    S = load(mat_file);
    data = S.(data_var);

    t_raw = data(:,1);
    Y = data(:,2:end);

    t_min = min(t_raw);
    t_max = max(t_raw);
    t = 2*(t_raw - t_min)/(t_max - t_min) - 1;

    n_samples = length(t);
    n_components = size(Y,2);

    best_degrees = zeros(n_components,1);
    coefficients = NaN(n_components,max_degree+1);
    component_rmse = zeros(n_components,1);
    Y_fit = zeros(size(Y));
    cv_rmse_by_component = NaN(n_components,max_degree);

    cv = cvpartition(n_samples,'KFold',K);

    for c = 1:n_components

        y = Y(:,c);
        degree_rmse = NaN(max_degree,1);
        ill_conditioned = false(max_degree,1);

        for deg = 1:max_degree

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

            degree_rmse(deg) = mean(fold_rmse,'omitnan');
            ill_conditioned(deg) = degree_bad;
        end

        valid = ~ill_conditioned;

        if any(valid)
            temp = degree_rmse;
            temp(~valid) = Inf;
            [~, best_deg] = min(temp);
        else
            error("No stable degree found for Orbit=%s, Component=%s", Orbit, component_names(c));
        end

        best_degrees(c) = best_deg;
        cv_rmse_by_component(c,:) = degree_rmse.';

        V_final = buildVandermonde(t,best_deg);
        a_final = V_final \ y;
        y_fit = V_final * a_final;

        Y_fit(:,c) = y_fit;
        coefficients(c,1:length(a_final)) = a_final';
        component_rmse(c) = sqrt(mean((y - y_fit).^2));
    end

    pos_error = Y(:,1:3) - Y_fit(:,1:3);
    vel_error = Y(:,4:6) - Y_fit(:,4:6);
    acc_error = Y(:,7:9) - Y_fit(:,7:9);

    position_rmse = sqrt(mean(sum(pos_error.^2,2)));
    velocity_rmse = sqrt(mean(sum(vel_error.^2,2)));
    acceleration_rmse = sqrt(mean(sum(acc_error.^2,2)));

    fprintf('\n=====================================================\n');
    fprintf('LSA RESULTS FOR %s\n', Orbit);
    fprintf('=====================================================\n');

    for c = 1:n_components
        fprintf('%3s | Best Degree = %2d | RMSE = %.6e\n', ...
            component_names(c), best_degrees(c), component_rmse(c));
    end

    fprintf('-----------------------------------------------------\n');
    fprintf('Position RMSE     = %.6e\n', position_rmse);
    fprintf('Velocity RMSE     = %.6e\n', velocity_rmse);
    fprintf('Acceleration RMSE = %.6e\n', acceleration_rmse);
    fprintf('=====================================================\n');

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

    nexttile;
    hold on;
    grid on;

    for c = 1:n_components
        plot(1:max_degree, cv_rmse_by_component(c,:), ...
            'LineWidth', 1.2, ...
            'DisplayName', component_names(c));
    end

    for c = 1:n_components
        bd = best_degrees(c);
        plot(bd, cv_rmse_by_component(c,bd), 'o', ...
            'MarkerSize', 6, ...
            'LineWidth', 1.2, ...
            'HandleVisibility','off');
    end

    xlabel('Polynomial Degree');
    ylabel('Val RMSE');
    title(sprintf('%s Orbit', Orbit));

    if o == 1
        legend('Location','eastoutside');
    end
end

set(gcf,'color','w');
sgtitle('Best-Fit Polynomial Degree Selection using 5-Fold Cross-Validation');

saveas(gcf, "LSA/All_Orbits_CV_RMSE_vs_Degree.png");

clearvars -except Final_Table

function V = buildVandermonde(t,degree)

n = length(t);
V = zeros(n,degree+1);

for k = 0:degree
    V(:,k+1) = t.^k;
end

end