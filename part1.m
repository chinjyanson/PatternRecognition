function part1(varargin)
    % Main function that calls others
    arg_length = length(varargin);
    showFigs = varargin{arg_length};
    for i = 1:arg_length-1
        displayEndEffectorPositions(varargin{i}, showFigs);
    end
end

function displayEndEffectorPositions(data, showFigs)
    if ~showFigs
        return;  % Exit function early if not showing figures
    end
    pos = data.end_effector_poses;

    figure
    plot3(pos(:,1), pos(:,2), pos(:,3))
    grid on
    axis equal

    xlabel('X')
    ylabel('Y')
    zlabel('Z')
    title('End-Effector Position Trajectory')

    figure
    subplot(3,1,1)
    plot(pos(:,4))
    ylabel('Roll')

    subplot(3,1,2)
    plot(pos(:,5))
    ylabel('Pitch')

    subplot(3,1,3)
    plot(pos(:,6))
    ylabel('Yaw')
    xlabel('Time (s)')
end