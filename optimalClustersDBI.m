function optimalK = optimalClustersDBI(data, maxK)
    % data: An NxM matrix with N data points and M features.
    % maxK: The maximum number of clusters to evaluate.
    
    % Initialize an array to store the DBI for each k
    dbIndices = zeros(1, maxK-1);

    for k = 2:maxK
        labels = kmeans(data, k);
        dbIndices(k-1) = daviesBouldinIndex(data, labels);
    end
    
    % Plotting DBI values
    figure;
    plot(2:maxK, dbIndices, 'o-', 'LineWidth', 2);
    xlabel('Number of Clusters (k)');
    ylabel('Davies-Bouldin Index');
    title('DBI for Different k Values');
    grid on;

    % Find the optimal number of clusters
    [~, optimalK] = min(dbIndices);
    optimalK = optimalK + 1;  % Adjusting for the indexing offset since we started at k=2
end

% DBI calculation function
function dbi = daviesBouldinIndex(data, labels)
    k = max(labels);
    clusterCentroids = arrayfun(@(i) mean(data(labels == i, :), 1), 1:k, 'UniformOutput', false);
    clusterCentroids = cell2mat(clusterCentroids');

    S = zeros(k, 1);
    for i = 1:k
        S(i) = mean(sqrt(sum(bsxfun(@minus, data(labels == i, :), clusterCentroids(i, :)).^2, 2)));
    end

    M = zeros(k, k);
    for i = 1:k
        for j = 1:k
            M(i, j) = sqrt(sum((clusterCentroids(i, :) - clusterCentroids(j, :)).^2));
        end
    end

    R = zeros(k, k);
    for i = 1:k
        for j = 1:k
            if i ~= j
                R(i, j) = (S(i) + S(j)) / M(i, j);
            end
        end
    end

    dbi = mean(max(R, [], 2));
end
