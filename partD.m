function partD(extractedData, showFigs)
    % partD  Linear Discriminant Analysis on extracted contact displacement data.
    % Usage: partD(extractedData, showFigs)
    %   extractedData - struct from partA containing contact data for all materials
    %   showFigs - boolean to control figure display

    if nargin < 1
        tmp = load('extractedData.mat');
        f = fieldnames(tmp);
        extractedData = tmp.(f{1});
    end
    if nargin < 2
        showFigs = true;
    end

    % Find oblong_TPU and oblong_rubber in extracted data
    fields = fieldnames(extractedData);
    data1 = [];
    data2 = [];

    for i = 1:numel(fields)
        d = extractedData.(fields{i});
        if ~isfield(d, 'name')
            continue;
        end
        n = lower(string(d.name));
        if contains(n, 'oblong') && contains(n, 'tpu')
            data1 = d;
        elseif contains(n, 'oblong') && contains(n, 'rubber')
            data2 = d;
        end
    end

    if isempty(data1) || isempty(data2)
        error('partD requires oblong_TPU and oblong_rubber data in extractedData.');
    end

    % Extract central papillae displacement (columns 10:12 for middle sensor)
    [X1, X2] = extractCentralDisplacement(data1, data2);

    % b. Visualize 3D displacement
    if showFigs
        visualize3DDisplacement(X1, X2, data1.name, data2.name);
    end

    % c. Apply LDA to all 2D combinations
    if showFigs
        applyLDA2DCombinations(X1, X2, data1.name, data2.name);
    end

    % d. Apply LDA to 3D displacement data
    applyLDA3D(X1, X2, data1.name, data2.name, showFigs);
end

function [X1, X2] = extractCentralDisplacement(data1, data2)
    % Central papillae (P4) displacement is columns 13:15 (D_X, D_Y, D_Z)
    % Layout: P0(1:3), P1(4:6), P2(7:9), P3(10:12), P4(13:15), P5(16:18), P6(19:21), P7(22:24), P8(25:27)
    % Supports both raw data (sensor_matrices_displacement) and extracted data (displacement)

    % Extract from data1
    if isfield(data1, 'displacement')
        X1 = data1.displacement(:, 13:15);
    elseif isfield(data1, 'sensor_matrices_displacement')
        X1 = data1.sensor_matrices_displacement(:, 13:15);
    else
        error('Input data1 is missing displacement data field.');
    end

    % Extract from data2
    if isfield(data2, 'displacement')
        X2 = data2.displacement(:, 13:15);
    elseif isfield(data2, 'sensor_matrices_displacement')
        X2 = data2.sensor_matrices_displacement(:, 13:15);
    else
        error('Input data2 is missing displacement data field.');
    end
end

function visualize3DDisplacement(X1, X2, name1, name2)
    figure('Name', 'Tactile Displacement 3D Scatter');
    hold on
    Utilities.plotByMaterial(X1(:,1), X1(:,2), X1(:,3), name1, ...
        'LineStyle', 'none', 'Marker', '.', 'MarkerSize', 15, 'DisplayName', char(name1));
    Utilities.plotByMaterial(X2(:,1), X2(:,2), X2(:,3), name2, ...
        'LineStyle', 'none', 'Marker', '.', 'MarkerSize', 15, 'DisplayName', char(name2));
    hold off
    grid on; axis equal
    xlabel('D_X'); ylabel('D_Y'); zlabel('D_Z');
    title('Central Papillae Tactile Displacement (3D)')
    legend('show', 'Location', 'bestoutside')
    view(3)
    
    fprintf('\nObservation (3D displacement):\n');
    fprintf('%s: mean [%.4f, %.4f, %.4f], std [%.4f, %.4f, %.4f]\n', ...
        name1, mean(X1), std(X1));
    fprintf('%s: mean [%.4f, %.4f, %.4f], std [%.4f, %.4f, %.4f]\n', ...
        name2, mean(X2), std(X2));
    fprintf('The two materials show %s in their displacement patterns.\n\n', ...
        determineOverlap(X1, X2));
end

function overlap = determineOverlap(X1, X2)
    mu1 = mean(X1, 1);
    mu2 = mean(X2, 1);
    dist = norm(mu1 - mu2);
    avgStd = mean([std(X1(:)), std(X2(:))]);
    
    if dist > 2 * avgStd
        overlap = "clear separation";
    elseif dist > avgStd
        overlap = "moderate separation";
    else
        overlap = "significant overlap";
    end
end

function applyLDA2DCombinations(X1, X2, name1, name2)
    % Apply LDA to all 2D combinations: (D_X, D_Y), (D_X, D_Z), (D_Y, D_Z)
    combos = {[1,2], [1,3], [2,3]};
    labels = ["D_X vs D_Y", "D_X vs D_Z", "D_Y vs D_Z"];
    axisLabels = {{'D_X', 'D_Y'}, {'D_X', 'D_Z'}, {'D_Y', 'D_Z'}};
    
    for c = 1:numel(combos)
        idx = combos{c};
        X1c = X1(:, idx);
        X2c = X2(:, idx);
        
        [w, threshold] = computeLDA2D(X1c, X2c);
        
        figure('Name', sprintf('LDA 2D | %s', labels(c)));
        hold on
        
        % Plot data points
        Utilities.plotByMaterial(X1c(:,1), X1c(:,2), name1, ...
            'LineStyle', 'none', 'Marker', '.', 'MarkerSize', 15, 'DisplayName', char(name1));
        Utilities.plotByMaterial(X2c(:,1), X2c(:,2), name2, ...
            'LineStyle', 'none', 'Marker', '.', 'MarkerSize', 15, 'DisplayName', char(name2));
        
        % Plot decision boundary
        plotDecisionBoundary2D(w, threshold);
        
        hold off
        grid on; axis equal
        xlabel(axisLabels{c}{1}); ylabel(axisLabels{c}{2});
        title(sprintf('LDA: %s', labels(c)))
        legend('show', 'Location', 'bestoutside')
        
        % Compute accuracy
        acc = computeAccuracy2D(X1c, X2c, w, threshold);
        fprintf('LDA %s: accuracy %.2f%%\n', labels(c), acc * 100);
    end
end

function [w, threshold] = computeLDA2D(X1, X2)
    % Compute LDA projection vector and threshold
    mu1 = mean(X1, 1)';
    mu2 = mean(X2, 1)';
    
    Sw = cov(X1, 1) + cov(X2, 1);  % Within-class scatter
    
    % Add small regularization for numerical stability
    Sw = Sw + eye(size(Sw)) * 1e-6;
    
    % LDA weight vector
    w = Sw \ (mu1 - mu2);
    w = w / norm(w);  % Normalize
    
    % Threshold: midpoint of projected means
    proj1 = X1 * w;
    proj2 = X2 * w;
    threshold = (mean(proj1) + mean(proj2)) / 2;
end

function plotDecisionBoundary2D(w, threshold)
    % Plot decision boundary line: w' * x = threshold
    % Rearrange to: w(2)*y = threshold - w(1)*x
    xlims = xlim;
    ylims = ylim;
    
    if abs(w(2)) > 1e-6
        x_line = linspace(xlims(1), xlims(2), 100);
        y_line = (threshold - w(1) * x_line) / w(2);
        plot(x_line, y_line, 'k-', 'LineWidth', 2, 'DisplayName', 'Decision boundary');
    else
        % Vertical line
        x_val = threshold / w(1);
        plot([x_val, x_val], ylims, 'k-', 'LineWidth', 2, 'DisplayName', 'Decision boundary');
    end
end

function acc = computeAccuracy2D(X1, X2, w, threshold)
    proj1 = X1 * w;
    proj2 = X2 * w;
    
    correct1 = sum(proj1 > threshold);
    correct2 = sum(proj2 <= threshold);
    
    acc = (correct1 + correct2) / (size(X1,1) + size(X2,1));
end

function applyLDA3D(X1, X2, name1, name2, showFigs)
    % d.i. Reduce to 2D and plot with LD and discrimination lines
    [w, threshold] = computeLDA3D(X1, X2);
    
    % Project to 2D using LDA direction and orthogonal direction
    w1 = w / norm(w);
    
    % Find orthogonal direction with maximum variance
    X_all = [X1; X2];
    X_centered = X_all - mean(X_all, 1);
    X_proj = X_centered - (X_centered * w1) * w1';
    [~, ~, V] = svd(X_proj, 'econ');
    w2 = V(:, 1);  % First PC of residuals
    
    % Project both datasets
    Y1 = [X1 * w1, X1 * w2];
    Y2 = [X2 * w1, X2 * w2];
    
    if showFigs
        figure('Name', 'LDA 3D reduced to 2D');
        hold on
        Utilities.plotByMaterial(Y1(:,1), Y1(:,2), name1, ...
            'LineStyle', 'none', 'Marker', '.', 'MarkerSize', 15, 'DisplayName', char(name1));
        Utilities.plotByMaterial(Y2(:,1), Y2(:,2), name2, ...
            'LineStyle', 'none', 'Marker', '.', 'MarkerSize', 15, 'DisplayName', char(name2));
        
        % Decision line (perpendicular to LD1 axis)
        ylims = ylim;
        plot([threshold, threshold], ylims, 'k-', 'LineWidth', 2, 'DisplayName', 'Decision line');
        
        hold off
        grid on; axis equal
        xlabel('LD1 (discriminant direction)'); ylabel('LD2 (orthogonal)');
        title('LDA: 3D to 2D projection with decision line')
        legend('show', 'Location', 'bestoutside')
    end
    
    % Compute accuracy
    proj1 = X1 * w1;
    proj2 = X2 * w1;
    correct1 = sum(proj1 > threshold);
    correct2 = sum(proj2 <= threshold);
    acc = (correct1 + correct2) / (size(X1,1) + size(X2,1));
    fprintf('\nLDA 3D accuracy: %.2f%%\n', acc * 100);
    
    % d.ii. Show 3D plot with discrimination plane
    if showFigs
        plot3DWithDiscriminationPlane(X1, X2, name1, name2, w, threshold);
    end
end

function [w, threshold] = computeLDA3D(X1, X2)
    mu1 = mean(X1, 1)';
    mu2 = mean(X2, 1)';
    
    Sw = cov(X1, 1) + cov(X2, 1);
    Sw = Sw + eye(size(Sw)) * 1e-6;
    
    w = Sw \ (mu1 - mu2);
    w = w / norm(w);
    
    proj1 = X1 * w;
    proj2 = X2 * w;
    threshold = (mean(proj1) + mean(proj2)) / 2;
end

function plot3DWithDiscriminationPlane(X1, X2, name1, name2, w, threshold)
    figure('Name', 'LDA 3D with discrimination plane');
    hold on

    Utilities.plotByMaterial(X1(:,1), X1(:,2), X1(:,3), name1, ...
        'LineStyle', 'none', 'Marker', '.', 'MarkerSize', 15, 'DisplayName', char(name1));
    Utilities.plotByMaterial(X2(:,1), X2(:,2), X2(:,3), name2, ...
        'LineStyle', 'none', 'Marker', '.', 'MarkerSize', 15, 'DisplayName', char(name2));

    % Plot discrimination plane: w' * x = threshold
    plotDiscriminationPlane(w, threshold, X1, X2);

    % Fit axes to data range with padding
    X_all = [X1; X2];
    margin = 0.55;
    for ax = 1:3
        rng_ax = max(X_all(:,ax)) - min(X_all(:,ax));
        pad = margin * rng_ax;
        switch ax
            case 1, xlim([min(X_all(:,1))-pad, max(X_all(:,1))+pad]);
            case 2, ylim([min(X_all(:,2))-pad, max(X_all(:,2))+pad]);
            case 3, zlim([min(X_all(:,3))-pad, max(X_all(:,3))+pad]);
        end
    end

    hold off
    grid on
    xlabel('D_X'); ylabel('D_Y'); zlabel('D_Z');
    title('LDA: 3D with discrimination plane')
    legend('show', 'Location', 'bestoutside')
    view(3)
end

function plotDiscriminationPlane(w, threshold, X1, X2)
    % Create plane: w(1)*x + w(2)*y + w(3)*z = threshold
    % Clip plane to data bounds so it doesn't extend beyond the box
    X_all = [X1; X2];
    x_range = [min(X_all(:,1)), max(X_all(:,1))];
    y_range = [min(X_all(:,2)), max(X_all(:,2))];
    z_range = [min(X_all(:,3)), max(X_all(:,3))];

    [X_grid, Y_grid] = meshgrid(linspace(x_range(1), x_range(2), 20), ...
                                 linspace(y_range(1), y_range(2), 20));

    if abs(w(3)) > 1e-6
        Z_grid = (threshold - w(1) * X_grid - w(2) * Y_grid) / w(3);
        % Clip plane to data z-range with padding
        z_pad = 0.5 * (z_range(2) - z_range(1));
        Z_grid(Z_grid < z_range(1) - z_pad) = NaN;
        Z_grid(Z_grid > z_range(2) + z_pad) = NaN;
        surf(X_grid, Y_grid, Z_grid, 'FaceAlpha', 0.3, 'EdgeColor', 'none', ...
            'FaceColor', 'k', 'DisplayName', 'Discrimination plane');
    else
        % Degenerate case: plot as line
        warning('Discrimination plane is degenerate (perpendicular to z-axis).');
    end
end