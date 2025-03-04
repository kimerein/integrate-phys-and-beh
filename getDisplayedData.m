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

% Process Scatter Objects
% Find scatter objects in the current axes.
% (Scatter objects in MATLAB are of type 'Scatter'.)
scatterObjs = findobj(ax, 'Type', 'Scatter');

% Initialize a cell array to store visible data points for each scatter object
visibleScatterData = cell(length(scatterObjs), 1);

% Loop through each scatter object
for k = 1:length(scatterObjs)
    % Get the data for this scatter object
    xData = scatterObjs(k).XData;
    yData = scatterObjs(k).YData;
    
    % Determine which points are within the current axis limits
    inView = (xData >= xLimits(1)) & (xData <= xLimits(2)) & ...
             (yData >= yLimits(1)) & (yData <= yLimits(2));
    
    % Extract only the visible data points
    visibleX = xData(inView);
    visibleY = yData(inView);
    
    % Retrieve the color data for the scatter object
    visibleColor = scatterObjs(k).MarkerEdgeColor;
    
    % Store the visible data points along with their color in the cell array
    visibleScatterData{k} = struct('x', visibleX, 'y', visibleY, 'color', visibleColor);
    
    % Optionally, display the number of visible points for this scatter object
    fprintf('Scatter %d: %d points visible\n', k, sum(inView));
end

% takeColors={[0 0.7500 0],[0 1 1],[0.8000 0.8000 0.8000],[1 0 0]};
% takeThisPointX=[]; takeThisPointY=[]; for i=1:length(visibleScatterData)
% isMatch = any(cellfun(@(c) isequal(c, visibleScatterData{i}.color), takeColors));
% if isMatch==true
% takeThisPointX=[takeThisPointX visibleScatterData{i}.x];
% takeThisPointY=[takeThisPointY visibleScatterData{i}.y];
% end
% end
% figure(); scatter(takeThisPointX-2.27604,takeThisPointY);

end
