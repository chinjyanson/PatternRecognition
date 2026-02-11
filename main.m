% Config
base = "PR_CW_mat/";
showFigs = true;
set(0, 'DefaultTextInterpreter', 'none');
set(0, 'DefaultLegendInterpreter', 'none');
set(0, 'DefaultAxesTickLabelInterpreter', 'none');

% Data Loading
cylinder_single = load(base + "cylinder_papillarray_single.mat");
cylinder_single.name = "cylinder_single";
cylinder_rubber = load(base + "cylinder_rubber_papillarray_single.mat");
cylinder_rubber.name = "cylinder_rubber";
cylinder_TPU    = load(base + "cylinder_TPU_papillarray_single.mat");
cylinder_TPU.name = "cylinder_TPU";
hexagon_single  = load(base + "hexagon_papillarray_single.mat");
hexagon_single.name = "hexagon_single";
hexagon_rubber  = load(base + "hexagon_rubber_papillarray_single.mat");
hexagon_rubber.name = "hexagon_rubber";
hexagon_TPU     = load(base + "hexagon_TPU_papillarray_single.mat");
hexagon_TPU.name = "hexagon_TPU";
oblong_single   = load(base + "oblong_papillarray_single.mat");
oblong_single.name = "oblong_single";
oblong_rubber   = load(base + "oblong_rubber_papillarray_single.mat");
oblong_rubber.name = "oblong_rubber";
oblong_TPU      = load(base + "oblong_TPU_papillarray_single.mat");
oblong_TPU.name = "oblong_TPU";
disp("Data loading completed")

% Part A - Visualization and contact data extraction
% Returns extracted contact data (force/displacement at peak contact moments)
extractedData = partA(cylinder_single, cylinder_rubber, cylinder_TPU, ...
    hexagon_single, hexagon_rubber, hexagon_TPU, ...
    oblong_single, oblong_rubber, oblong_TPU, showFigs);

% Part B - PCA on extracted force data
partB(extractedData, showFigs)

% Part C - t-SNE on extracted force data
partC(extractedData, showFigs)

% Part D - LDA on extracted displacement data
partD(extractedData, showFigs)

% Part E - Clustering on extracted force/displacement data
partE(extractedData, showFigs)

% Part F - GMM on extracted displacement data
partF(extractedData, showFigs)

% Part G - Bagging on extracted displacement data
%partG(extractedData, showFigs)
