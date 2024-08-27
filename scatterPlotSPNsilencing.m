function scatterPlotSPNsilencing(row2Array, row4Array)
    figure();
    % Check if the first argument has two rows
    if size(row2Array, 1) ~= 2
        error('The first input array must have exactly two rows.');
    end
    
    % Check if the second argument has four rows
    if size(row4Array, 1) ~= 4
        error('The second input array must have exactly four rows.');
    end

    % Find the indices where the third row of the four-row array is 1
    validIndices = (row4Array(3, :) == 1);

    % Extract the corresponding columns from the two-row array
    xData = row2Array(1, validIndices);
    yData = row2Array(2, validIndices);

    % Create the scatter plot
    scatter(xData, yData, 'filled');
    
    % Label the axes
    xlabel('Top row of first argument');
    ylabel('Second row of first argument');
    
    % Set a title for the plot
    title('Scatter plot of filtered data');
    
    % Display the grid for better visualization
    grid on;
end
