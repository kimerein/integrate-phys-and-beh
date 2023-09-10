function optimalK = silhouetteAnalysis(data, maxK)
    % data: An NxM matrix with N data points and M features.
    % maxK: The maximum number of clusters to evaluate.

    % Remove rows with NaN values
    data(any(isnan(data), 2), :) = [];
    
    avgSilhouette = zeros(1, maxK-1); % Average silhouette values

    for k = 2:maxK % Silhouette is not defined for k=1
        clusterLabels = kmeans(data, k);
        silhouetteValues = silhouette(data, clusterLabels);
        avgSilhouette(k-1) = mean(silhouetteValues);
    end

    % Plotting
    figure;
    plot(2:maxK, avgSilhouette, 'o-', 'LineWidth', 2);
    xlabel('Number of Clusters (k)');
    ylabel('Average Silhouette Value');
    title('Silhouette Analysis for Optimal k');
    grid on;

    % Return the k with maximum average silhouette value
    [~, optimalK] = max(avgSilhouette);
    optimalK = optimalK + 1; % Adjusting for the indexing offset due to starting at k=2
end

