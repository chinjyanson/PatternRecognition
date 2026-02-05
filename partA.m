function extractedData = partA(varargin)
    % partA  Visualize end-effector trajectories and force/torque data for all object shapes.
    %        Also extracts contact data at maximum force moments for use in later parts.
    % Usage: extractedData = partA(cyl_normal, cyl_rubber, cyl_tpu, hex_normal, hex_rubber, hex_tpu,
    %                              oblong_normal, oblong_rubber, oblong_tpu, showFigs)
    % Returns: extractedData - struct containing force and displacement data at contact peaks

    arg_length = length(varargin);
    showFigs = varargin{arg_length};

    % Parse data and separate into shapes
    cylinders = {};
    hexagons = {};
    oblongs = {};

    for i = 1:arg_length-1
        data = varargin{i};
        if isfield(data, 'name')
            n = lower(string(data.name));
            if contains(n, 'cylinder')
                cylinders{end+1} = data; %#ok<AGROW>
            elseif contains(n, 'hexagon')
                hexagons{end+1} = data; %#ok<AGROW>
            elseif contains(n, 'oblong')
                oblongs{end+1} = data; %#ok<AGROW>
            end
        end
    end

    % A.1: Display end-effector positions for cylinder and hexagon PLA (normal) only
    for i = 1:numel(cylinders)
        data = cylinders{i};
        if contains(lower(string(data.name)), 'normal')
            displayEndEffectorPositions(data, data.name, showFigs);
        end
    end
    for i = 1:numel(hexagons)
        data = hexagons{i};
        if contains(lower(string(data.name)), 'normal')
            displayEndEffectorPositions(data, data.name, showFigs);
        end
    end

    % A.2: Force/torque peaks and contact extraction for all objects
    allData = [cylinders, hexagons, oblongs];
    extractedData = struct();

    for i = 1:numel(allData)
        data = allData{i};
        materialName = data.name;
        ft_peaks = obtainForceTorquePeaks(data);
        displayForceTorque(data, ft_peaks, materialName, showFigs);

        % Extract contact data at force peaks for future parts (A.2.c)
        contactData = extractContactData(data, ft_peaks);
        fieldName = matlab.lang.makeValidName(materialName);
        extractedData.(fieldName) = contactData;
    end

    % A.3: 3D scatter plots of EXTRACTED contact force data from middle papillae
    % One plot per shape, with different colors for each material
    if showFigs
        plotExtractedForceByShape(cylinders, extractedData, 'Cylinder');
        plotExtractedForceByShape(hexagons, extractedData, 'Hexagon');
        plotExtractedForceByShape(oblongs, extractedData, 'Oblong');
    end

    disp('Part A completed. Contact data extracted for future parts.');
end

function plotExtractedForceByShape(dataList, extractedData, shapeLabel)
    % A.3.a: Plot 3D scatter of middle papillae force for extracted contacts
    % P4 (middle papillae) columns: 13:15
    if isempty(dataList)
        return;
    end

    figure('Name', sprintf('%s - Middle Papillae Force (Extracted Contacts)', shapeLabel));
    hold on

    for i = 1:numel(dataList)
        data = dataList{i};
        fieldName = matlab.lang.makeValidName(data.name);
        if isfield(extractedData, fieldName)
            contactData = extractedData.(fieldName);
            F = contactData.force(:, 13:15);  % P4 columns
            Utilities.plotByMaterial(F(:,1), F(:,2), F(:,3), data.name, ...
                'LineStyle', 'none', 'Marker', '.', 'MarkerSize', 15, ...
                'DisplayName', char(data.name));
        end
    end

    hold off
    grid on
    axis equal
    xlabel('F_X')
    ylabel('F_Y')
    zlabel('F_Z')
    title(sprintf('%s Objects - Middle Papillae (P4) Force at Contact Peaks', shapeLabel))
    legend('show', 'Location', 'bestoutside')
end

function contactData = extractContactData(data, ft_peaks)
    % Extract sensor data at maximum force contact moments
    % Uses Fz peaks (column 3) as the primary indicator of contact
    % Extracts ONE sample per contact (at the peak moment)

    contactData = struct();
    contactData.name = data.name;

    % Use Fz (column 3) peaks as contact indicators - this is the normal force
    fz_peak_indices = ft_peaks{3}.indices;

    if isempty(fz_peak_indices)
        warning('No Fz peaks found for %s. Using all data.', data.name);
        contactData.force = data.sensor_matrices_force;
        contactData.displacement = data.sensor_matrices_displacement;
        contactData.ft_values = data.ft_values;
        contactData.peak_indices = [];
        return;
    end

    % Extract data at each peak index (one row per contact)
    contactData.force = data.sensor_matrices_force(fz_peak_indices, :);
    contactData.displacement = data.sensor_matrices_displacement(fz_peak_indices, :);
    contactData.ft_values = data.ft_values(fz_peak_indices, :);
    contactData.peak_indices = fz_peak_indices;

    fprintf('  %s: Extracted %d contact samples (one per peak)\n', ...
        data.name, length(fz_peak_indices));
end

function displayEndEffectorPositions(data, material, showFigs)
    if ~showFigs
        return;  % Exit function early if not showing figures
    end
    pos = data.end_effector_poses;

    figure
    Utilities.plotByMaterial(pos(:,1), pos(:,2), pos(:,3), material)
    grid on
    axis equal

    xlabel('X')
    ylabel('Y')
    zlabel('Z')
    title('End-Effector Position Trajectory')

    figure
    subplot(3,1,1)
    Utilities.plotByMaterial(1:length(pos(:,4)), pos(:,4), material)
    ylabel('Roll')

    subplot(3,1,2)
    Utilities.plotByMaterial(1:length(pos(:,5)), pos(:,5), material)
    ylabel('Pitch')

    subplot(3,1,3)
    Utilities.plotByMaterial(1:length(pos(:,6)), pos(:,6), material)
    ylabel('Yaw')
    xlabel('Time (s)')
end

function displayForceTorque(data, peaks, material, showFigs)
    if ~showFigs
        return;
    end

    ft_data = data.ft_values;

    % Force plots in first figure
    figure;
    F = ft_data(:,1:3);

    for j = 1:3
        subplot(3,1,j)
        Utilities.plotByMaterial(1:length(F(:,j)), F(:,j), material)
        hold on
        if ~isempty(peaks{j}.indices)
            plot(peaks{j}.indices, peaks{j}.values, 'r^', 'MarkerSize', 8, 'DisplayName', 'Peaks')
        end
        hold off
        ylabel(sprintf('F%s (N)', char(119+j)))  % Fx, Fy, Fz
        legend
    end
    sgtitle(sprintf('Force Values - %s', material))

    % Torque plots in second figure
    figure;
    T = ft_data(:,4:6);

    for j = 1:3
        subplot(3,1,j)
        Utilities.plotByMaterial(1:length(T(:,j)), T(:,j), material)
        hold on
        if ~isempty(peaks{j+3}.indices)
            plot(peaks{j+3}.indices, peaks{j+3}.values, 'r^', 'MarkerSize', 8, 'DisplayName', 'Peaks')
        end
        hold off
        ylabel(sprintf('T%s (Nm)', char(119+j)))  % Tx, Ty, Tz
        legend
    end
    sgtitle(sprintf('Torque Values - %s', material))
end

function peaks = obtainForceTorquePeaks(data)
    ft_data = data.ft_values;
    
    % Initialize cell arrays to store peaks and troughs
    peaks = cell(6, 1);
    
    labels = {'Fx', 'Fy', 'Fz', 'Tx', 'Ty', 'Tz'};
    
    for i = 1:6
        signal = ft_data(:, i);
        
        % Find peaks (local maxima)
        [peak_values, peak_indices] = findLocalMaxima(signal);
        peaks{i}.values = peak_values;
        peaks{i}.indices = peak_indices;
        
        fprintf('%s: %d peaks\n', labels{i}, length(peak_values));
    end
end

function [values, indices] = findLocalMaxima(signal)
    % Find local maxima (peaks) in signal across 30 data points
    n = length(signal);
    indices = [];
    values = [];
    window = 1000;  % Check ±15 points (30 total)
    
    for i = window+1:n-window
        % Check if current point is maximum in the window
        is_peak = true;
        for j = -window:window
            if j ~= 0 && signal(i) <= signal(i+j)
                is_peak = false;
                break;
            end
        end
        
        if is_peak
            indices = [indices; i];
            values = [values; signal(i)];
        end
    end
end

