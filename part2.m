function part2(varargin)
    % part2  PCA analysis on force data for cylinder and oblong objects.
    % Usage: part2(cyl_normal, cyl_rubber, cyl_tpu, oblong_normal,
    %              oblong_rubber, oblong_tpu, showFigs)

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
end

function runGroupPCA(dataList, groupLabel, showFigs)
    if isempty(dataList)
        return;
    end

    [X, materials] = collectForceData(dataList);
    mu = mean(X, 1);
    sigma = std(X, 0, 1);
    sigma(sigma == 0) = eps;  % avoid divide-by-zero
    Xz = (X - mu) ./ sigma;  % manual z-score (no toolbox)

    % PCA without Statistics Toolbox: eigendecomposition of covariance
    C = cov(Xz, 1);  % normalize by N (like MATLAB pca 'Rows','pairwise' default with standardization)
    [V, D] = eig(C);
    latent = max(diag(D), 0);                % guard tiny negatives
    [latent, idx] = sort(latent, 'descend'); % sort descending
    coeff = V(:, idx);
    score = Xz * coeff;

    if ~showFigs
        return;
    end

    plotStandardized3D(Xz, materials, coeff, latent, groupLabel);
    plotScores2D(score, materials, groupLabel);
    plotComponentNumberLines(score, materials, groupLabel);
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
