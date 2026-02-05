function partG(extractedData, showFigs)
    % partG  Bagging (Bootstrap Aggregation) on displacement data with PCA.
    % Usage: partG(extractedData, showFigs)
    %   extractedData - struct from partA containing contact data for all materials
    %   showFigs - boolean to control figure display
    %
    % Applies bagging to displacement data from all nine papillae,
    % previously processed with PCA

    if nargin < 2
        showFigs = true;
    end

    fprintf('\n=== Part G: Bagging (Bootstrap Aggregation) ===\n');

    % Collect displacement data from ALL 9 papillae for all objects
    fields = fieldnames(extractedData);
    allData = {};

    for i = 1:numel(fields)
        d = extractedData.(fields{i});
        if isfield(d, 'name')
            allData{end+1} = d; %#ok<AGROW>
        end
    end

    [X, labels] = collectAllPapillaeDisplacement(allData);

    fprintf('Total samples: %d\n', size(X,1));
    fprintf('Features (27 displacement values from 9 papillae): %d\n', size(X,2));

    % Apply PCA for dimensionality reduction (as mentioned in spec: "previously processed with PCA")
    fprintf('\nApplying PCA for dimensionality reduction...\n');
    [X_pca, ~, nComponents] = applyPCA(X);
    fprintf('Reduced to %d principal components (95%% variance explained)\n', nComponents);

    % Split data into training and test sets (80/20)
    rng(42);  % For reproducibility
    cv = cvpartition(labels, 'HoldOut', 0.2);
    X_train = X_pca(cv.training, :);
    X_test = X_pca(cv.test, :);
    y_train = labels(cv.training);
    y_test = labels(cv.test);

    fprintf('\nTraining samples: %d\n', length(y_train));
    fprintf('Test samples: %d\n', length(y_test));

    % G.1.a: Specify number of bags/trees
    nTrees = 50;  % Number of decision trees
    fprintf('\nNumber of trees (bags): %d\n', nTrees);
    fprintf('Reason: 50 trees provides good balance between accuracy and computation.\n');
    fprintf('        More trees reduce variance but with diminishing returns after ~50-100.\n');

    % Train bagged ensemble of decision trees (Random Forest)
    fprintf('\nTraining bagged decision tree ensemble...\n');
    baggedModel = TreeBagger(nTrees, X_train, y_train, ...
        'Method', 'classification', ...
        'OOBPrediction', 'on', ...
        'OOBPredictorImportance', 'on');

    % G.1.b: Visualize two decision trees
    if showFigs
        visualizeDecisionTrees(baggedModel);
    end

    % G.1.c: Run model on test data and display confusion matrix
    y_pred = predict(baggedModel, X_test);

    % Convert to categorical for confusionmat
    y_test_cat = categorical(y_test);
    y_pred_cat = categorical(y_pred);

    % Get unique classes in consistent order
    classes = categories(y_test_cat);

    if showFigs
        plotConfusionMatrix(y_test_cat, y_pred_cat, classes);
    end

    % Calculate and report accuracy
    accuracy = sum(y_pred_cat == y_test_cat) / length(y_test_cat);
    fprintf('\nOverall Test Accuracy: %.2f%%\n', accuracy * 100);

    % Per-class accuracy
    fprintf('\nPer-class accuracy:\n');
    for i = 1:length(classes)
        classIdx = y_test_cat == classes{i};
        classAcc = sum(y_pred_cat(classIdx) == y_test_cat(classIdx)) / sum(classIdx);
        fprintf('  %s: %.2f%%\n', classes{i}, classAcc * 100);
    end

    % Out-of-bag error
    oobErr = oobError(baggedModel);
    fprintf('\nOut-of-bag error: %.2f%%\n', oobErr(end) * 100);

    % G.1.d: Discussion points
    fprintf('\n=== Discussion ===\n');
    fprintf('G.1.d Analysis:\n');
    fprintf('  - Misclassifications often occur between materials with similar physical properties\n');
    fprintf('  - Rubber and TPU may be confused as both are soft/compliant materials\n');
    fprintf('  - PCA helps by reducing noise and focusing on principal variations\n');
    fprintf('  - PCA may also lose some discriminative information in lower components\n');

    % Plot OOB error vs number of trees
    if showFigs
        plotOOBError(baggedModel);
        plotFeatureImportance(baggedModel, nComponents);
    end

    fprintf('\nPart G completed.\n');
end

function [X, labels] = collectAllPapillaeDisplacement(dataList)
    % Collect displacement data from ALL 9 papillae (27 features)
    X = [];
    labels = strings(0,1);

    for i = 1:numel(dataList)
        d = dataList{i};
        if isfield(d, 'displacement')
            D = d.displacement(:, 1:27);  % All 9 papillae
        elseif isfield(d, 'sensor_matrices_displacement')
            D = d.sensor_matrices_displacement(:, 1:27);
        else
            error('Input data is missing displacement field.');
        end
        X = [X; D]; %#ok<AGROW>
        labels = [labels; repmat(string(d.name), size(D,1), 1)]; %#ok<AGROW>
    end
end

function [X_pca, explained, nComponents] = applyPCA(X)
    % Apply PCA and retain components explaining 95% variance
    [~, score, ~, ~, explained] = pca(X);

    % Find number of components for 95% variance
    cumVar = cumsum(explained);
    nComponents = find(cumVar >= 95, 1);
    if isempty(nComponents)
        nComponents = size(score, 2);
    end

    X_pca = score(:, 1:nComponents);
end

function visualizeDecisionTrees(baggedModel)
    % G.1.b: Visualize two decision trees from the ensemble
    figure('Name', 'Decision Tree 1');
    view(baggedModel.Trees{1}, 'Mode', 'graph');
    title('Decision Tree 1 from Bagged Ensemble')

    figure('Name', 'Decision Tree 2');
    view(baggedModel.Trees{2}, 'Mode', 'graph');
    title('Decision Tree 2 from Bagged Ensemble')
end

function plotConfusionMatrix(y_true, y_pred, classes)
    % G.1.c: Plot confusion matrix
    figure('Name', 'Confusion Matrix');

    % Compute confusion matrix
    C = confusionmat(y_true, y_pred);

    % Plot as heatmap
    heatmap(classes, classes, C, ...
        'Title', 'Confusion Matrix', ...
        'XLabel', 'Predicted Class', ...
        'YLabel', 'True Class', ...
        'ColorbarVisible', 'on');

    % Also display as text
    fprintf('\nConfusion Matrix:\n');
    fprintf('%-20s', 'True\\Predicted');
    for i = 1:length(classes)
        fprintf('%-15s', classes{i});
    end
    fprintf('\n');

    for i = 1:length(classes)
        fprintf('%-20s', classes{i});
        for j = 1:length(classes)
            fprintf('%-15d', C(i,j));
        end
        fprintf('\n');
    end
end

function plotOOBError(baggedModel)
    % Plot out-of-bag error vs number of trees
    figure('Name', 'OOB Error vs Trees');

    oobErr = oobError(baggedModel);
    plot(1:length(oobErr), oobErr, 'b-', 'LineWidth', 2);

    xlabel('Number of Trees');
    ylabel('Out-of-Bag Classification Error');
    title('OOB Error vs Number of Trees');
    grid on

    fprintf('\nFinal OOB error (with all trees): %.4f\n', oobErr(end));
end

function plotFeatureImportance(baggedModel, ~)
    % Plot feature importance from the bagged model
    figure('Name', 'Feature Importance');

    importance = baggedModel.OOBPermutedPredictorDeltaError;

    bar(importance);
    xlabel('Principal Component');
    ylabel('Importance (OOB Permuted Delta Error)');
    title('Feature Importance by Principal Component');
    grid on

    % Label top features
    [~, sortIdx] = sort(importance, 'descend');
    fprintf('\nTop 5 most important principal components:\n');
    for i = 1:min(5, length(sortIdx))
        fprintf('  PC%d: %.4f\n', sortIdx(i), importance(sortIdx(i)));
    end
end
