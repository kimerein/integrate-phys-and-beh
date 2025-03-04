function getDisplayedData()

% Get the current axes and its limits
ax = gca;
xLimits = ax.XLim;
yLimits = ax.YLim;

% Find all line objects in the current axes
lineObjs = findobj(ax, 'Type', 'line');

% Initialize a cell array to store visible data points for each line
visibleData = cell(length(lineObjs), 1);

% Loop through each line object
for k = 1:length(lineObjs)
    % Get the data for this line
    xData = lineObjs(k).XData;
    yData = lineObjs(k).YData;
    
    if yLimits(1)==0
        yLimits(1)=0-0.0001;
    end

    % Determine which points are within the current axis limits
    inView = (xData >= xLimits(1)) & (xData <= xLimits(2)) & ...
             (yData >= yLimits(1)) & (yData <= yLimits(2));
    
    % Extract only the visible data points
    visibleX = xData(inView);
    visibleY = yData(inView);

    % Store the visible data points in the cell array as a structure
    visibleData{k} = struct('x', visibleX, 'y', visibleY);
    
    % Display number of visible points for this line
    fprintf('Line %d: %d points visible\n', k, sum(inView));
end

end
