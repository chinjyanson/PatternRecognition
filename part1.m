function part1(varargin)
    % Main function that calls others
    arg_length = length(varargin);
    showFigs = varargin{arg_length};
    for i = 1:arg_length-1
        materialName = inputname(i);  % Get variable name from caller
        displayEndEffectorPositions(varargin{i}, materialName, showFigs);
        [ft_peaks, ft_troughs] = obtainForceTorquePeaks(varargin{i});
        displayForceTorque(varargin{i}, ft_peaks, ft_troughs, materialName, showFigs);
        display3DForceDataMiddlePapillae(varargin{i}, materialName, showFigs);
    end
end

function display3DForceDataMiddlePapillae(data, material, showFigs)
    if ~showFigs
        return;  % Exit function early if not showing figures
    end
    forceData = data.sensor_matrices_force;

    figure
    Utilities.plotByMaterial(forceData(:,10), forceData(:,11), forceData(:,12), material)
    grid on
    axis equal

    xlabel('X')
    ylabel('Y')
    zlabel('Z')
    title('Force Data of Middle Papillae')
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

function displayForceTorque(data, peaks, troughs, material, showFigs)
    if ~showFigs
        return;
    end

    ft_data = data.ft_values;  % Extract the force/torque data from the structure
    
    % Force plots in first figure
    figure;
    F = ft_data(:,1:3);
    
    for j = 1:3
        subplot(3,1,j)
        Utilities.plotByMaterial(1:length(F(:,j)), F(:,j), material)
        hold on
        % Plot peaks
        if ~isempty(peaks{j}.indices)
            plot(peaks{j}.indices, peaks{j}.values, 'r^', 'MarkerSize', 8, 'DisplayName', 'Peaks')
        end
        % Plot troughs
        if ~isempty(troughs{j}.indices)
            plot(troughs{j}.indices, troughs{j}.values, 'gv', 'MarkerSize', 8, 'DisplayName', 'Troughs')
        end
        hold off
        ylabel(sprintf('F%s (N)', char(119+j)))  % Fx, Fy, Fz
        legend
    end
    sgtitle('Force Values')

    % Torque plots in second figure
    figure;
    T = ft_data(:,4:6);

    for j = 1:3
        subplot(3,1,j)
        Utilities.plotByMaterial(1:length(T(:,j)), T(:,j), material)
        hold on
        % Plot peaks
        if ~isempty(peaks{j+3}.indices)
            plot(peaks{j+3}.indices, peaks{j+3}.values, 'r^', 'MarkerSize', 8, 'DisplayName', 'Peaks')
        end
        % Plot troughs
        if ~isempty(troughs{j+3}.indices)
            plot(troughs{j+3}.indices, troughs{j+3}.values, 'gv', 'MarkerSize', 8, 'DisplayName', 'Troughs')
        end
        hold off
        ylabel(sprintf('T%s (Nm)', char(119+j)))  % Tx, Ty, Tz
        legend
    end
    sgtitle('Torque Values')
end

function [peaks, troughs] = obtainForceTorquePeaks(data)
    ft_data = data.ft_values;
    
    % Initialize cell arrays to store peaks and troughs
    peaks = cell(6, 1);
    troughs = cell(6, 1);
    
    labels = {'Fx', 'Fy', 'Fz', 'Tx', 'Ty', 'Tz'};
    
    for i = 1:6
        signal = ft_data(:, i);
        
        % Find peaks (local maxima)
        [peak_values, peak_indices] = findLocalMaxima(signal);
        peaks{i}.values = peak_values;
        peaks{i}.indices = peak_indices;
        
        % Find troughs (local minima)
        [trough_values, trough_indices] = findLocalMinima(signal);
        troughs{i}.values = trough_values;
        troughs{i}.indices = trough_indices;
        
        fprintf('%s: %d peaks, %d troughs\n', labels{i}, length(peak_values), length(trough_values));
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

function [values, indices] = findLocalMinima(signal)
    % Find local minima (troughs) in signal across 30 data points
    n = length(signal);
    indices = [];
    values = [];
    window = 1000;  % Check ±15 points (30 total)
    
    for i = window+1:n-window
        % Check if current point is minimum in the window
        is_trough = true;
        for j = -window:window
            if j ~= 0 && signal(i) >= signal(i+j)
                is_trough = false;
                break;
            end
        end
        
        if is_trough
            indices = [indices; i];
            values = [values; signal(i)];
        end
    end
end
