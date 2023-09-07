function clusterIdx=clusterAfterDTW(data,numClusters,maxsamp)

numClusters=numClusters+1;
wereNan=all(isnan(data),2);

% Assuming you have a matrix 'data' where each row is a time series

% Step 1: Compute pairwise DTW distances
numSeries = size(data, 1);
dtwDistances = zeros(numSeries);

for i = 1:numSeries
    for j = i+1:numSeries
        dtwDistances(i, j) = dtw(data(i, :), data(j, :), maxsamp);
        dtwDistances(j, i) = dtwDistances(i, j);  % The matrix is symmetric
    end
end

dtwDistances(isnan(dtwDistances))=0;

% Step 2: Apply kmedoids using the DTW distances
clusterIdx = kmeans(dtwDistances, numClusters);

clusterIdx(wereNan)=nan;

% 'clusterIdx' contains the cluster index for each time series
% 'medoid' contains the index of the representative time series for each cluster

end