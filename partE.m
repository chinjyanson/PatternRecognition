function partE(extractedData, showFigs)
    % partE  Clustering analysis on extracted contact data.
    % Usage: partE(extractedData, showFigs)
    %   extractedData - struct from partA containing contact data for all materials
    %   showFigs - boolean to control figure display
    %
    % Uses oblong objects with displacement data from central papillae (P4)

    if nargin < 2
        showFigs = true;
    end

    % Get oblong data (can also use hexagon by changing this)
    fields = fieldnames(extractedData);
    oblongs = {};

    for i = 1:numel(fields)
        d = extractedData.(fields{i});
        if ~isfield(d, 'name')
            continue;
        end
        n = lower(string(d.name));
        if contains(n, 'oblong')
            oblongs{end+1} = d; %#ok<AGROW>
        end
    end

    if isempty(oblongs)
        error('No oblong data found in extractedData.');
    end

    % Collect displacement data from central papillae (P4) for all 3 materials
    % P4 columns: 13:15
    [X, materials] = collectDisplacementData(oblongs);

    fprintf('\n=== Part E: Clustering ===\n');
    fprintf('Using displacement data from central papillae (P4)\n');
    fprintf('Total samples: %d\n\n', size(X,1));

    % E.1.a: Scatter plot with ground truth colors
    if showFigs
        plotGroundTruth(X, materials, 'Oblong');
    end

    % E.1.b: Apply K-means clustering with Euclidean distance
    k = 3;  % 3 materials = 3 clusters
    [clusterIdx1, centroids1] = applyKMeans(X, k, 'sqeuclidean');

    if showFigs
        plotClusteringResult(X, materials, clusterIdx1, centroids1, ...
            'K-Means (Euclidean)', 'Oblong');
    end

    % Analyze clustering vs ground truth
    analyzeClusteringAccuracy(materials, clusterIdx1, 'K-Means (Euclidean)');

    % E.1.c: Change distance metric and repeat
    % Using cityblock (Manhattan) distance
    [clusterIdx2, centroids2] = applyKMeans(X, k, 'cityblock');

    if showFigs
        plotClusteringResult(X, materials, clusterIdx2, centroids2, ...
            'K-Means (Manhattan)', 'Oblong');
    end

    analyzeClusteringAccuracy(materials, clusterIdx2, 'K-Means (Manhattan)');

    fprintf('\nPart E completed.\n');
end

function [X, materials] = collectDisplacementData(dataList)
    % Collect displacement data from central papillae (P4) - columns 13:15
    X = [];
    materials = strings(0,1);

    for i = 1:numel(dataList)
        d = dataList{i};
        if isfield(d, 'displacement')
            D = d.displacement(:, 13:15);
        elseif isfield(d, 'sensor_matrices_displacement')
            D = d.sensor_matrices_displacement(:, 13:15);
        else
            error('Input data is missing displacement field.');
        end
        X = [X; D]; %#ok<AGROW>
        materials = [materials; repmat(string(d.name), size(D,1), 1)]; %#ok<AGROW>
    end
end

function plotGroundTruth(X, materials, shapeLabel)
    % E.1.a: Plot scatter with ground truth material colors
    figure('Name', sprintf('%s - Ground Truth (Displacement)', shapeLabel));
    hold on

    mats = unique(materials);
    for i = 1:numel(mats)
        idx = materials == mats(i);
        Utilities.plotByMaterial(X(idx,1), X(idx,2), X(idx,3), mats(i), ...
            'LineStyle', 'none', 'Marker', '.', 'MarkerSize', 15, ...
            'DisplayName', char(mats(i)));
    end

    hold off
    grid on
    axis equal
    xlabel('D_X'); ylabel('D_Y'); zlabel('D_Z');
    title(sprintf('%s - Ground Truth Labels (Displacement)', shapeLabel))
    legend('show', 'Location', 'bestoutside')
end

function [clusterIdx, centroids] = applyKMeans(X, k, distanceMetric)
    % Apply K-means clustering
    fprintf('Applying K-means with %s distance...\n', distanceMetric);

    rng(42);  % For reproducibility
    [clusterIdx, centroids] = kmeans(X, k, ...
        'Distance', distanceMetric, ...
        'MaxIter', 1000, ...
        'Replicates', 10);

    fprintf('  Clusters found: %d\n', k);
end

function plotClusteringResult(X, materials, clusterIdx, centroids, methodName, shapeLabel)
    % Plot clustering result with cluster assignments shown by marker shape
    figure('Name', sprintf('%s - %s', shapeLabel, methodName));
    hold on

    markers = {'o', 's', 'd', '^', 'v', '>', '<', 'p', 'h'};
    colors = lines(max(clusterIdx));
    mats = unique(materials);

    % Plot each point with color = material, shape = cluster
    for m = 1:numel(mats)
        for c = 1:max(clusterIdx)
            idx = (materials == mats(m)) & (clusterIdx == c);
            if any(idx)
                % Use material color from Utilities, marker from cluster
                scatter3(X(idx,1), X(idx,2), X(idx,3), 50, ...
                    'Marker', markers{c}, ...
                    'MarkerEdgeColor', colors(c,:), ...
                    'MarkerFaceColor', 'none', ...
                    'DisplayName', sprintf('%s (Cluster %d)', char(mats(m)), c));
            end
        end
    end

    % Plot centroids
    scatter3(centroids(:,1), centroids(:,2), centroids(:,3), 200, 'k', 'x', ...
        'LineWidth', 3, 'DisplayName', 'Centroids');

    hold off
    grid on
    axis equal
    xlabel('D_X'); ylabel('D_Y'); zlabel('D_Z');
    title(sprintf('%s - %s Clustering', shapeLabel, methodName))
    legend('show', 'Location', 'bestoutside')

    % Also create a 2D view for clarity
    figure('Name', sprintf('%s - %s (2D)', shapeLabel, methodName));

    % Plot with color = ground truth material
    subplot(1,2,1)
    hold on
    for i = 1:numel(mats)
        idx = materials == mats(i);
        Utilities.plotByMaterial(X(idx,1), X(idx,2), mats(i), ...
            'LineStyle', 'none', 'Marker', '.', 'MarkerSize', 15, ...
            'DisplayName', char(mats(i)));
    end
    hold off
    grid on; axis equal
    xlabel('D_X'); ylabel('D_Y');
    title('Ground Truth')
    legend('show', 'Location', 'best')

    % Plot with color = cluster assignment
    subplot(1,2,2)
    hold on
    for c = 1:max(clusterIdx)
        idx = clusterIdx == c;
        scatter(X(idx,1), X(idx,2), 30, colors(c,:), 'filled', ...
            'DisplayName', sprintf('Cluster %d', c));
    end
    scatter(centroids(:,1), centroids(:,2), 200, 'k', 'x', ...
        'LineWidth', 3, 'DisplayName', 'Centroids');
    hold off
    grid on; axis equal
    xlabel('D_X'); ylabel('D_Y');
    title(methodName)
    legend('show', 'Location', 'best')

    sgtitle(sprintf('%s Objects - Clustering Comparison', shapeLabel))
end

function analyzeClusteringAccuracy(materials, clusterIdx, methodName)
    % Analyze how well clusters correspond to ground truth materials
    mats = unique(materials);
    k = max(clusterIdx);

    fprintf('\n%s - Cluster composition:\n', methodName);

    % Create confusion-like matrix
    for c = 1:k
        fprintf('  Cluster %d: ', c);
        clusterMask = (clusterIdx == c);
        total = sum(clusterMask);

        for m = 1:numel(mats)
            matMask = (materials == mats(m));
            count = sum(clusterMask & matMask);
            pct = 100 * count / total;
            fprintf('%s: %d (%.1f%%)  ', mats(m), count, pct);
        end
        fprintf('\n');
    end

    % Calculate purity (best-case assignment accuracy)
    purity = 0;
    for c = 1:k
        clusterMask = (clusterIdx == c);
        maxCount = 0;
        for m = 1:numel(mats)
            matMask = (materials == mats(m));
            count = sum(clusterMask & matMask);
            if count > maxCount
                maxCount = count;
            end
        end
        purity = purity + maxCount;
    end
    purity = purity / length(clusterIdx);
    fprintf('  Purity: %.2f%%\n', purity * 100);
end
