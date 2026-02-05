function part2(varargin)
    % part2  PCA analysis on force data for cylinder and oblong objects.
    % Usage: part2(cyl_normal, cyl_rubber, cyl_tpu, hex_normal, hex_rubber, hex_tpu,
    %              oblong_normal, oblong_rubber, oblong_tpu, showFigs)
    % Note: hexagon data is accepted but only cylinder and oblong are analyzed.

    if nargin == 0
        error('part2 expects data structures as input.');
    end

    if islogical(varargin{end})
        showFigs = varargin{end};
        dataArgs = varargin(1:end-1);
    else
        showFigs = true;
        dataArgs = varargin;
    end

    % Split inputs into cylinder and oblong groups by name
    cylinders = {};
    oblongs = {};
    for i = 1:numel(dataArgs)
        d = dataArgs{i};
        if ~isfield(d, 'name')
            error('Each input must have a .name field indicating material.');
        end
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
    runAllPapillaePCA(dataArgs, showFigs);
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
    % Middle sensor force is stored in sensor_matrices_force columns 10:12
    X = [];
    materials = strings(0,1);
    for i = 1:numel(dataList)
        d = dataList{i};
        if ~isfield(d, 'sensor_matrices_force')
            error('Input data is missing sensor_matrices_force field.');
        end
        F = d.sensor_matrices_force(:, 10:12);  % Fx,Fy,Fz of middle sensor
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
    X = [];
    materials = strings(0,1);
    for i = 1:numel(dataList)
        d = dataList{i};
        if ~isfield(d, 'sensor_matrices_force')
            error('Input data is missing sensor_matrices_force field.');
        end
        Ffull = d.sensor_matrices_force;
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
