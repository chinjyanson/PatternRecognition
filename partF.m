function partF(extractedData, showFigs)
    % partF  Gaussian Mixture Model analysis on displacement data.
    % Usage: partF(extractedData, showFigs)
    %   extractedData - struct from partA containing contact data for all materials
    %   showFigs - boolean to control figure display
    %
    % Classifies objects into material types using displacement data

    if nargin < 2
        showFigs = true;
    end

    % Select an object shape (using oblong here)
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

    % Collect displacement data from central papillae (P4) - columns 13:15
    [X, materials] = collectDisplacementData(oblongs);

    % Use D_X and D_Z (columns 1 and 3) for better material separation
    X2D = X(:, [1, 3]);

    fprintf('\n=== Part F: Gaussian Mixture Model ===\n');
    fprintf('Using 2D displacement data (D_X, D_Z) from central papillae\n');
    fprintf('Total samples: %d\n\n', size(X2D,1));

    % F.1.a: 2D scatter plot with ground truth colors
    if showFigs
        plotGroundTruth2D(X2D, materials, 'Oblong');
    end

    % F.1.b: Fit GMM with 3 components and plot contour
    nComponents = 3;
    gmmModel = fitGMM(X2D, nComponents);

    if showFigs
        plotGMMContour(X2D, materials, gmmModel, 'Oblong');
    end

    % F.1.c: Plot GMM as 3D surface
    if showFigs
        plotGMM3DSurface(X2D, gmmModel, 'Oblong');
    end

    % F.1.d: Assign hard clusters and replot
    clusterIdx = cluster(gmmModel, X2D);

    if showFigs
        plotHardClusters(X2D, materials, clusterIdx, 'Oblong');
    end

    % F.1.e: Compare clusters to ground truth
    analyzeGMMClusters(materials, clusterIdx);

    fprintf('\nPart F completed.\n');
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

function plotGroundTruth2D(X, materials, shapeLabel)
    % F.1.a: 2D scatter plot with ground truth material colors
    figure('Name', sprintf('%s - GMM Ground Truth', shapeLabel));
    hold on

    mats = unique(materials);
    for i = 1:numel(mats)
        idx = materials == mats(i);
        Utilities.plotByMaterial(X(idx,1), X(idx,2), mats(i), ...
            'LineStyle', 'none', 'Marker', '.', 'MarkerSize', 15, ...
            'DisplayName', char(mats(i)));
    end

    hold off
    grid on
    axis equal
    xlabel('D_X (Displacement) (m)');
    ylabel('D_Z (Displacement) (m)');
    title(sprintf('%s - Ground Truth Material Labels', shapeLabel))
    legend('show', 'Location', 'bestoutside')
end

function gmmModel = fitGMM(X, nComponents)
    % F.1.b: Fit Gaussian Mixture Model
    fprintf('Fitting GMM with %d components...\n', nComponents);

    rng(42);  % For reproducibility
    options = statset('MaxIter', 1000);

    % Use minimal regularization to capture true data covariance
    % Lower regularization = tighter, more accurate contours
    gmmModel = fitgmdist(X, nComponents, ...
        'CovarianceType', 'full', ...
        'SharedCovariance', false, ...
        'RegularizationValue', 1e-6, ...
        'Options', options, ...
        'Replicates', 5);  % Try multiple starts for better fit

    fprintf('  Converged: %s\n', mat2str(gmmModel.Converged));
    fprintf('  NegLogLikelihood: %.4f\n', gmmModel.NegativeLogLikelihood);

    % Print component means and proportions
    fprintf('  Component means:\n');
    for k = 1:nComponents
        fprintf('    C%d: [%.4f, %.4f], weight=%.2f%%\n', k, ...
            gmmModel.mu(k,1), gmmModel.mu(k,2), gmmModel.ComponentProportion(k)*100);
    end
end

function plotGMMContour(X, materials, gmmModel, shapeLabel)
    % F.1.b: Plot GMM probability density contours with scatter overlay
    % Shows separate contour regions for each Gaussian component
    figure('Name', sprintf('%s - GMM Contour Plot', shapeLabel));

    % Create grid for contour evaluation - use data range with small margin
    dataRangeX = max(X(:,1)) - min(X(:,1));
    dataRangeY = max(X(:,2)) - min(X(:,2));
    marginX = max(0.02, dataRangeX * 0.3);  % 30% margin or minimum 0.02
    marginY = max(0.02, dataRangeY * 0.3);

    xRange = [min(X(:,1)) - marginX, max(X(:,1)) + marginX];
    yRange = [min(X(:,2)) - marginY, max(X(:,2)) + marginY];

    [xGrid, yGrid] = meshgrid(linspace(xRange(1), xRange(2), 200), ...
                               linspace(yRange(1), yRange(2), 200));
    gridPoints = [xGrid(:), yGrid(:)];

    nComponents = gmmModel.NumComponents;

    % Compute PDF for each component separately (weighted by mixing proportion)
    componentPDFs = zeros(size(gridPoints, 1), nComponents);
    for k = 1:nComponents
        mu_k = gmmModel.mu(k, :);
        Sigma_k = gmmModel.Sigma(:, :, k);
        weight_k = gmmModel.ComponentProportion(k);

        % Compute multivariate normal PDF for this component
        diff = gridPoints - mu_k;
        mahalDist = sum((diff / Sigma_k) .* diff, 2);
        detSigma = det(Sigma_k);
        componentPDFs(:, k) = weight_k * (2*pi)^(-1) * detSigma^(-0.5) * exp(-0.5 * mahalDist);
    end

    % Compute total GMM PDF (sum of all components)
    totalPDF = sum(componentPDFs, 2);
    totalPDFGrid = reshape(totalPDF, size(xGrid));

    hold on

    % Draw unified contour lines for total GMM PDF
    contour(xGrid, yGrid, totalPDFGrid, 20, 'LineWidth', 1, ...
        'HandleVisibility', 'off');
    colormap(parula);
    cb = colorbar;
    ylabel(cb, 'Probability Density');

    % Mark component means with labels
    for k = 1:nComponents
        mu_k = gmmModel.mu(k, :);
        scatter(mu_k(1), mu_k(2), 200, 'k', 'p', 'filled', ...
            'MarkerEdgeColor', 'w', 'LineWidth', 2, 'HandleVisibility', 'off');
        text(mu_k(1) + 0.005, mu_k(2) + 0.01, sprintf('C%d', k), ...
            'FontSize', 12, 'FontWeight', 'bold', 'Color', 'k', ...
            'BackgroundColor', 'w');
    end

    % Overlay scatter plot with ground truth colors
    mats = unique(materials);
    for i = 1:numel(mats)
        idx = materials == mats(i);
        Utilities.plotByMaterial(X(idx,1), X(idx,2), mats(i), ...
            'LineStyle', 'none', 'Marker', '.', 'MarkerSize', 15, ...
            'DisplayName', char(mats(i)));
    end

    hold off
    grid on
    xlabel('D_X (Displacement) (m)');
    ylabel('D_Z (Displacement) (m)');
    title(sprintf('%s - F.b: GMM Probability Density Contours with Data Overlay', shapeLabel))
    legend('show', 'Location', 'bestoutside')
end

function plotGMM3DSurface(X, gmmModel, shapeLabel)
    % F.1.c: Plot GMM as 3D surface
    figure('Name', sprintf('%s - GMM 3D Surface', shapeLabel));

    % Create grid
    margin = 0.1;
    xRange = [min(X(:,1)) - margin, max(X(:,1)) + margin];
    yRange = [min(X(:,2)) - margin, max(X(:,2)) + margin];

    [xGrid, yGrid] = meshgrid(linspace(xRange(1), xRange(2), 80), ...
                               linspace(yRange(1), yRange(2), 80));
    gridPoints = [xGrid(:), yGrid(:)];

    % Evaluate GMM PDF
    pdfValues = pdf(gmmModel, gridPoints);
    pdfGrid = reshape(pdfValues, size(xGrid));

    % Plot surface
    surf(xGrid, yGrid, pdfGrid, 'EdgeColor', 'none', 'FaceAlpha', 0.8);
    colormap('jet');
    colorbar;

    hold on
    % Mark the means
    scatter3(gmmModel.mu(:,1), gmmModel.mu(:,2), ...
        pdf(gmmModel, gmmModel.mu), 200, 'k', 'x', 'LineWidth', 3);
    hold off

    xlabel('D_X (Displacement) (m)');
    ylabel('D_Z (Displacement) (m)');
    zlabel('Probability Density');
    title(sprintf('%s - GMM Probability Density Surface', shapeLabel))

    % Set isometric viewing angle
    view(45, 30);
    grid on
end

function plotHardClusters(X, materials, clusterIdx, shapeLabel)
    % F.1.d: Assign hard clusters and replot
    figure('Name', sprintf('%s - GMM Hard Clusters', shapeLabel));

    colors = lines(max(clusterIdx));

    % Plot side by side
    subplot(1,2,1)
    hold on
    mats = unique(materials);
    for i = 1:numel(mats)
        idx = materials == mats(i);
        Utilities.plotByMaterial(X(idx,1), X(idx,2), mats(i), ...
            'LineStyle', 'none', 'Marker', '.', 'MarkerSize', 15, ...
            'DisplayName', char(mats(i)));
    end
    hold off
    grid on; axis equal
    xlabel('D_X (m)'); ylabel('D_Z (m)');
    title('Ground Truth')
    legend('show', 'Location', 'best')

    subplot(1,2,2)
    hold on
    for c = 1:max(clusterIdx)
        idx = clusterIdx == c;
        scatter(X(idx,1), X(idx,2), 30, colors(c,:), 'filled', ...
            'DisplayName', sprintf('Cluster %d', c));
    end
    hold off
    grid on; axis equal
    xlabel('D_X (m)'); ylabel('D_Z (m)');
    title('GMM Hard Clusters')
    legend('show', 'Location', 'best')

    sgtitle(sprintf('%s - Ground Truth vs GMM Clusters', shapeLabel))
end

function analyzeGMMClusters(materials, clusterIdx)
    % F.1.e: Compare clusters to ground truth
    mats = unique(materials);
    k = max(clusterIdx);

    fprintf('\nGMM Cluster vs Ground Truth Analysis:\n');

    % Create confusion-like matrix
    fprintf('\nCluster composition:\n');
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

    % Calculate purity
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
    fprintf('\nClustering Purity: %.2f%%\n', purity * 100);

    fprintf('\nInterpretation:\n');
    fprintf('  - If clusters align well with materials, displacement alone can classify materials\n');
    fprintf('  - If clusters mix materials, shear displacement may not fully separate materials\n');
end
