function partB(extractedData, showFigs)
    % partB  PCA analysis on extracted contact force data for cylinder and oblong objects.
    % Usage: partB(extractedData, showFigs)
    %   extractedData - struct from partA containing contact data for all materials
    %   showFigs - boolean to control figure display
    % Note: hexagon data is accepted but only cylinder and oblong are analyzed.

    if nargin < 2
        showFigs = true;
    end

    % Get field names and separate by shape
    fields = fieldnames(extractedData);
    cylinders = {};
    oblongs = {};
    allData = {};

    for i = 1:numel(fields)
        d = extractedData.(fields{i});
        if ~isfield(d, 'name')
            continue;
        end
        allData{end+1} = d; %#ok<AGROW>
        n = lower(string(d.name));
        if contains(n, 'cylinder')
            cylinders{end+1} = d; %#ok<AGROW>
        elseif contains(n, 'oblong')
            oblongs{end+1} = d; %#ok<AGROW>
        end
    end

    runGroupPCA(cylinders, 'Cylinder', showFigs);
    runGroupPCA(oblongs, 'Oblong', showFigs);

    % Part 3: all nine papillae (all objects, one plot per shape)
    runAllPapillaePCA(allData, showFigs);
end

function summary = runGroupPCA(dataList, groupLabel, showFigs)
    if isempty(dataList)
        summary = struct();
        return;
    end

    [X, materials] = collectForceData(dataList);
    [Xz, coeff, score, latent] = computePCA(X);

    if ~showFigs
        summary = summarizePCA(score, materials, latent, groupLabel);
        return;
    end

    plotStandardized3D(Xz, materials, coeff, latent, groupLabel);
    plotScores2D(score, materials, groupLabel);
    plotComponentNumberLines(score, materials, groupLabel);
    summary = summarizePCA(score, materials, latent, groupLabel);
end

function [X, materials] = collectForceData(dataList)
    % Middle papillae (P4) force is stored in columns 13:15
    % Layout: P0(1:3), P1(4:6), P2(7:9), P3(10:12), P4(13:15), P5(16:18), P6(19:21), P7(22:24), P8(25:27)
    % Supports both raw data (sensor_matrices_force) and extracted data (force)
    X = [];
    materials = strings(0,1);
    for i = 1:numel(dataList)
        d = dataList{i};
        % Check for extracted contact data format first, then raw format
        if isfield(d, 'force')
            F = d.force(:, 13:15);  % Fx,Fy,Fz of middle sensor (P4)
        elseif isfield(d, 'sensor_matrices_force')
            F = d.sensor_matrices_force(:, 13:15);
        else
            error('Input data is missing force data field.');
        end
        X = [X; F]; %#ok<AGROW>
        materials = [materials; repmat(string(d.name), size(F,1), 1)]; %#ok<AGROW>
    end
end

function runAllPapillaePCA(dataArgs, showFigs)
    if isempty(dataArgs)
        return;
    end

    shapes = ["cylinder", "hexagon", "oblong"];
    labels = ["Cylinder", "Hexagon", "Oblong"];
    for s = 1:numel(shapes)
        list = {};
        for i = 1:numel(dataArgs)
            d = dataArgs{i};
            if isfield(d, 'name') && contains(lower(string(d.name)), shapes(s))
                list{end+1} = d; %#ok<AGROW>
            end
        end
        if ~isempty(list)
            runShapeAllPapillaePCA(list, labels(s), showFigs);
        end
    end
end

function runShapeAllPapillaePCA(dataList, groupLabel, showFigs)
    [X, materials] = collectAllPapillaeForceData(dataList);
    [~, ~, score, latent] = computePCA(X);

    if ~showFigs
        return;
    end

    figure('Name', groupLabel + " | All Papillae | PCA 2D scores");
    hold on
    plotMaterialScatter2D(score(:,1), score(:,2), materials);
    hold off
    grid on; axis equal
    xlabel('PC1 score'); ylabel('PC2 score');
    title(groupLabel + " | all papillae | 2D PCA projection")
    legend('show', 'Location', 'bestoutside')

    plotScree(latent, groupLabel + " | all papillae");
end

function [Xz, coeff, score, latent] = computePCA(X)
    mu = mean(X, 1);
    sigma = std(X, 0, 1);
    sigma(sigma == 0) = eps;  % avoid divide-by-zero
    Xz = (X - mu) ./ sigma;  % manual z-score standardization

    % Use Statistics Toolbox pca function
    [coeff, score, latent] = pca(Xz, 'Centered', false);
end

function summary = summarizePCA(score, materials, latent, groupLabel)
    explained = latent ./ max(sum(latent), eps);
    first2 = sum(explained(1:min(2, numel(explained))));

    mats = unique(materials);
    centroids = zeros(numel(mats), 2);
    for i = 1:numel(mats)
        idx = materials == mats(i);
        centroids(i,:) = mean(score(idx,1:2), 1);
    end

    if size(centroids,1) >= 2
        dsum = 0; cnt = 0;
        for i = 1:size(centroids,1)-1
            for j = i+1:size(centroids,1)
                dsum = dsum + norm(centroids(i,:) - centroids(j,:));
                cnt = cnt + 1;
            end
        end
        meanDist = dsum / cnt;
    else
        meanDist = NaN;
    end

    summary = struct(...
        'groupLabel', groupLabel, ...
        'explainedFirst2', first2, ...
        'meanCentroidDist2D', meanDist, ...
        'numSamples', size(score,1));
end

function [X, materials] = collectAllPapillaeForceData(dataList)
    % Supports both raw data (sensor_matrices_force) and extracted data (force)
    X = [];
    materials = strings(0,1);
    for i = 1:numel(dataList)
        d = dataList{i};
        % Check for extracted contact data format first, then raw format
        if isfield(d, 'force')
            Ffull = d.force;
        elseif isfield(d, 'sensor_matrices_force')
            Ffull = d.sensor_matrices_force;
        else
            error('Input data is missing force data field.');
        end
        if size(Ffull,2) >= 27
            F = Ffull(:, 1:27);
        else
            F = Ffull;
        end
        X = [X; F]; %#ok<AGROW>
        materials = [materials; repmat(string(d.name), size(F,1), 1)]; %#ok<AGROW>
    end
end

function plotScree(latent, titleLabel)
    figure('Name', titleLabel + " | Scree");
    explained = latent ./ max(sum(latent), eps);
    plot(1:numel(explained), 100 * explained, '-o', 'LineWidth', 1.5);
    grid on
    xlabel('Principal component');
    ylabel('Variance explained (%)');
    title(titleLabel + " | scree plot")
end

function plotStandardized3D(Xz, materials, coeff, latent, groupLabel)
    figure('Name', groupLabel + " | Standardized Force with PCs");
    hold on
    plotMaterialScatter3D(Xz(:,1), Xz(:,2), Xz(:,3), materials);

    % Overlay principal component directions
    scale = sqrt(latent(:));
    nComp = min(3, numel(scale));
    for i = 1:nComp
        quiver3(0, 0, 0, coeff(1,i)*scale(i), coeff(2,i)*scale(i), coeff(3,i)*scale(i), ...
            'k', 'LineWidth', 1.5, 'MaxHeadSize', 1, 'DisplayName', sprintf('PC%d', i));
    end
    hold off
    grid on; axis equal
    xlabel('Fx (z-scored)'); ylabel('Fy (z-scored)'); zlabel('Fz (z-scored)');
    title(groupLabel + " | Standardized force & principal components")
    legend('show', 'Location', 'bestoutside')
end

function plotScores2D(score, materials, groupLabel)
    figure('Name', groupLabel + " | PCA 2D scores");
    hold on
    plotMaterialScatter2D(score(:,1), score(:,2), materials);
    hold off
    grid on; axis equal
    xlabel('PC1 score'); ylabel('PC2 score');
    title(groupLabel + " | 2D PCA projection")
    legend('show', 'Location', 'bestoutside')
end

function plotComponentNumberLines(score, materials, groupLabel)
    nComp = size(score,2);
    figure('Name', groupLabel + " | PCA component number lines");
    for c = 1:nComp
        subplot(nComp,1,c)
        y = zeros(size(score,1),1) + c;  % position this component on its own line
        plotMaterialScatter2D(score(:,c), y, materials);
        grid on
        yticks(c)
        yticklabels("PC" + c)
        xlabel('Score value')
        if c == 1
            title(groupLabel + " | distribution across components (1D number lines)")
        end
    end
end

function plotMaterialScatter3D(x, y, z, materials)
    mats = unique(materials);
    for i = 1:numel(mats)
        idx = materials == mats(i);
        if any(idx)
            Utilities.plotByMaterial(x(idx), y(idx), z(idx), mats(i), ...
                'LineStyle', 'none', 'Marker', '.', 'MarkerSize', 8, 'DisplayName', char(mats(i)));
            hold on
        end
    end
end

function plotMaterialScatter2D(x, y, materials)
    mats = unique(materials);
    for i = 1:numel(mats)
        idx = materials == mats(i);
        if any(idx)
            Utilities.plotByMaterial(x(idx), y(idx), mats(i), ...
                'LineStyle', 'none', 'Marker', '.', 'MarkerSize', 8, 'DisplayName', char(mats(i)));
            hold on
        end
    end
end
