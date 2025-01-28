function [w]=plotNeuralProjection(matrix1, matrix2, times1, times2, time_window, passin_w)
% plotNeuralProjection Projects and plots neural activity projections for two conditions.
%
% This function takes two neural response matrices corresponding to different conditions,
% each with its own time vector, and a specified time window. It finds the optimal direction
% that best separates the two conditions by averaging the neural activity over the given time window.
% Each time point's activity is then projected onto this direction, and the projections
% are plotted over time for both conditions.
%
% Inputs:
%   matrix1    - N x T1 matrix for Condition 1 (N neurons x T1 timepoints)
%   matrix2    - N x T2 matrix for Condition 2 (N neurons x T2 timepoints)
%   times1     - 1 x T1 vector of time points corresponding to the columns of matrix1
%   times2     - 1 x T2 vector of time points corresponding to the columns of matrix2
%   time_window - 1 x 2 vector specifying the [start_time, end_time] for
%   averaging (or a cell array of 2 of these vectors, if want to use
%   different time windows for matrix1 and matrix2)
%
% Example usage:
%   plotNeuralProjection(cond1_matrix, cond2_matrix, time_vector1, time_vector2, [0.5, 1.5]);

    % Input Validation
    % Ensure both matrices have the same number of neurons
    if size(matrix1, 1) ~= size(matrix2, 1)
        error('Both matrices must have the same number of neurons (rows).');
    end
    
    % Ensure time vectors match the number of columns in their respective matrices
    if length(times1) ~= size(matrix1, 2)
        error('Length of times1 vector must match the number of columns in matrix1.');
    end
    if length(times2) ~= size(matrix2, 2)
        error('Length of times2 vector must match the number of columns in matrix2.');
    end
    
    time_window2=[];
    if iscell(time_window)
        if length(time_window)~=2
            error('time_window must be a two-element cell array (for matrix1 time window, then matrix2 time window), or a vector: [start_time, end_time]');
        else
            temp=time_window;
            time_window=temp{1};
            time_window2=temp{2};
        end
    else
        % Ensure time_window is a two-element vector
        if length(time_window) ~= 2
            error('time_window must be a two-element vector: [start_time, end_time].');
        end
    end
    
    % Ensure time_window is ordered correctly
    if time_window(1) > time_window(2)
        error('time_window should be in the format [start_time, end_time] with start_time <= end_time.');
    end

    % Identify Timepoints Within the Specified Window for Each Condition
    idx_window1 = times1 >= time_window(1) & times1 <= time_window(2);
    if ~isempty(time_window2)
        idx_window2 = times2 >= time_window3(1) & times2 <= time_window2(2);
    else
        idx_window2 = times2 >= time_window(1) & times2 <= time_window(2);
    end
    
    if ~any(idx_window1)
        warning('No timepoints found within the specified time window for Condition 1.');
    end
    if ~any(idx_window2)
        warning('No timepoints found within the specified time window for Condition 2.');
    end
    
    % Combine indices to handle cases where one condition has no points in the window
    if any(idx_window1) && any(idx_window2)
        % Average neural activity over the time window for each condition
        mean1 = mean(matrix1(:, idx_window1), 2, 'omitnan'); % N x 1 vector for Condition 1
        mean2 = mean(matrix2(:, idx_window2), 2, 'omitnan'); % N x 1 vector for Condition 2
    elseif any(idx_window1)
        % Only Condition 1 has points in the window
        mean1 = mean(matrix1(:, idx_window1), 2, 'omitnan');
        mean2 = mean(matrix2, 2, 'omitnan'); % Use overall mean for Condition 2
        warning('Only Condition 1 has timepoints within the specified window. Using overall mean for Condition 2.');
    elseif any(idx_window2)
        % Only Condition 2 has points in the window
        mean1 = mean(matrix1, 2, 'omitnan'); % Use overall mean for Condition 1
        mean2 = mean(matrix2(:, idx_window2), 2, 'omitnan');
        warning('Only Condition 2 has timepoints within the specified window. Using overall mean for Condition 1.');
    else
        % Neither condition has points in the window
        error('No timepoints found within the specified time window for either condition.');
    end
    
    % Determine the Separation Direction
    % The optimal direction is the difference between the mean vectors
    w = mean1 - mean2;

    % Zero any nans in w and before projecting
    w(isnan(w))=0;
    matrix1(isnan(matrix1))=0;
    matrix2(isnan(matrix2))=0;
    
    % Check if the separation direction has non-zero magnitude
    w_norm = norm(w);
    if w_norm == 0
        error('The separation direction has zero magnitude. Check your data and time window.');
    end
    
    % Normalize the Separation Direction to Unit Length
    w = w / w_norm;
    
    % If passed in w, use this w instead of comparison
    if ~isempty(passin_w)
        w=passin_w;
    end

    % Project Each Timepoint's Activity Onto the Separation Direction
    % For Condition 1
    proj1 = w' * matrix1; % 1 x T1 vector of projections
    % For Condition 2
    proj2 = w' * matrix2; % 1 x T2 vector of projections
    
    % Plotting the Projection Magnitudes Over Time
    figure;
    hold on;
    
    % Plot Condition 1 Projections
    plot(times1, proj1, 'r', 'LineWidth', 1.5, 'DisplayName', 'Condition 1');
    
    % Plot Condition 2 Projections
    plot(times2, proj2, 'b', 'LineWidth', 1.5, 'DisplayName', 'Condition 2');
    
    % Highlight the Averaging Time Window
    yLimits = ylim;
    plot([time_window(1), time_window(1)], yLimits, 'k--', 'LineWidth', 1);
    plot([time_window(2), time_window(2)], yLimits, 'k--', 'LineWidth', 1);
    legend('Location', 'best');
    
    % Enhance Plot Aesthetics
    xlabel('Time');
    ylabel('Projection Magnitude');
    title('Neural Activity Projection onto Separation Direction');
    grid on;
    hold off;
    
    % Additional Output (Optional)
    % You can return the separation vector and projections if needed
    % For example:
    % varargout = {w, proj1, proj2};
end