% Config
base = "PR_CW_mat/";
showFigs = true;

% Data Loading
cylinder_normal = load(base + "cylinder_papillarray_single.mat");
cylinder_rubber = load(base + "cylinder_rubber_papillarray_single.mat");
cylinder_TPU    = load(base + "cylinder_TPU_papillarray_single.mat");
hexagon_normal  = load(base + "hexagon_papillarray_single.mat");
hexagon_rubber  = load(base + "hexagon_rubber_papillarray_single.mat");
hexagon_TPU     = load(base + "hexagon_TPU_papillarray_single.mat");
oblong_normal   = load(base + "oblong_papillarray_single.mat");
oblong_rubber   = load(base + "oblong_rubber_papillarray_single.mat");
oblong_TPU      = load(base + "oblong_TPU_papillarray_single.mat");
disp("Data loading completed")

% Part 1 
part1(cylinder_normal, cylinder_rubber, cylinder_TPU, ...
    hexagon_normal, hexagon_rubber, hexagon_TPU, showFigs)

% Part 2


