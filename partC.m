function partC(extractedData, showFigs)
	% partC  t-SNE on extracted contact force data.
	% Usage: partC(extractedData, showFigs)
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

	if exist('tsne', 'file') ~= 2
		error('t-SNE requires the Statistics and Machine Learning Toolbox (tsne function not found).');
	end

	% Get field names and separate by shape
	fields = fieldnames(extractedData);
	cylinders = {};
	hexagons = {};
	oblongs = {};

	for i = 1:numel(fields)
		d = extractedData.(fields{i});
		if ~isfield(d, 'name')
			continue;
		end
		n = lower(string(d.name));
		if contains(n, 'cylinder')
			cylinders{end+1} = d; %#ok<AGROW>
		elseif contains(n, 'hexagon')
			hexagons{end+1} = d; %#ok<AGROW>
		elseif contains(n, 'oblong')
			oblongs{end+1} = d; %#ok<AGROW>
		end
	end

	perplexities = [5, 15];

	runTSNEGroup(cylinders, 'Cylinder', perplexities, showFigs);
	% Repeat for either hexagon or oblong (choose oblong here)
	runTSNEGroup(oblongs, 'Oblong', perplexities, showFigs);
	% runTSNEGroup(hexagons, 'Hexagon', perplexities, showFigs);
end

function summary = runTSNEGroup(dataList, groupLabel, perplexities, showFigs)
	summary = struct();
	if isempty(dataList)
		return;
	end

	[X, materials] = collectForceData(dataList);
	[Xz] = standardizeFeatures(X);

	losses = zeros(numel(perplexities), 1);
	meanDist = zeros(numel(perplexities), 1);

	for p = 1:numel(perplexities)
		rng(0);
		[Y, loss] = tsne(Xz, 'NumDimensions', 2, 'Perplexity', perplexities(p), ...
			'NumPCAComponents', min(50, size(Xz,2)), 'Verbose', 0);
		losses(p) = loss;
		meanDist(p) = meanCentroidDistance(Y, materials);

		if showFigs
			figure('Name', groupLabel + " | t-SNE (perplexity=" + perplexities(p) + ")");
			hold on
			plotMaterialScatter2D(Y(:,1), Y(:,2), materials);
			hold off
			grid on; axis equal
			xlabel('t-SNE 1'); ylabel('t-SNE 2');
			title(groupLabel + " | t-SNE (perplexity=" + perplexities(p) + ")")
			legend('show', 'Location', 'bestoutside')
		end
	end

	summary = struct(...
		'groupLabel', groupLabel, ...
		'perplexities', perplexities(:), ...
		'losses', losses, ...
		'meanCentroidDist2D', meanDist, ...
		'numSamples', size(Xz,1));

	fprintf('\n%s t-SNE losses:\n', groupLabel);
	for p = 1:numel(perplexities)
		fprintf('  Perplexity %d -> loss %.4f\n', perplexities(p), losses(p));
	end
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

function [Xz] = standardizeFeatures(X)
	mu = mean(X, 1);
	sigma = std(X, 0, 1);
	sigma(sigma == 0) = eps;
	Xz = (X - mu) ./ sigma;
end

function d = meanCentroidDistance(Y, materials)
	mats = unique(materials);
	centroids = zeros(numel(mats), 2);
	for i = 1:numel(mats)
		idx = materials == mats(i);
		centroids(i,:) = mean(Y(idx,1:2), 1);
	end
	if size(centroids,1) < 2
		d = NaN;
		return;
	end
	dsum = 0; cnt = 0;
	for i = 1:size(centroids,1)-1
		for j = i+1:size(centroids,1)
			dsum = dsum + norm(centroids(i,:) - centroids(j,:));
			cnt = cnt + 1;
		end
	end
	d = dsum / cnt;
end

function plotMaterialScatter2D(x, y, materials)
	mats = unique(materials);
	for i = 1:numel(mats)
		idx = materials == mats(i);
		if any(idx)
			Utilities.plotByMaterial(x(idx), y(idx), mats(i), ...
				'LineStyle', 'none', 'Marker', '.', 'MarkerSize', 15, 'DisplayName', char(mats(i)));
			hold on
		end
	end
end
