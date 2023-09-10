function optimalK = optimalClustersGapStatistic(data, maxK, B)
    % data: NxM data matrix
    % maxK: Maximum number of clusters to consider
    % B: Number of reference datasets for bootstrap sampling

    rng('default');  % For reproducibility

    % Remove rows containing NaN values
    data = data(~any(isnan(data), 2), :);

    % Compute the sum of squared distances for the actual dataset
    actualDists = arrayfun(@(k) sum(kmeans(data, k, 'EmptyAction', 'singleton', 'Replicates', 10).^2), 1:maxK);

    % Generate B reference datasets and compute their distances
    refDists = zeros(B, maxK);
    for i = 1:B
        refData = arrayfun(@(z) min(data(:, z)) + range(data(:, z)) * rand(size(data, 1), 1), 1:size(data, 2), 'UniformOutput', false);
        refData = horzcat(refData{:});
        for k = 1:maxK
            refDists(i, k) = sum(kmeans(refData, k, 'EmptyAction', 'singleton', 'Replicates', 10).^2);
        end
    end

    % Compute the Gap statistic
    gaps = mean(log(refDists) - log(repmat(actualDists, B, 1)));
    errs = std(log(refDists), 0, 1) * sqrt(1 + 1/B);

    % Find the optimal number of clusters
    optimalK = find(gaps > [0, gaps(1:end-1)] + [0, errs(1:end-1)], 1, 'last');

    % Plot the Gap values
    figure;
    bar(1:maxK, gaps);
    hold on;
    errorbar(1:maxK, gaps, errs, 'k', 'LineStyle', 'none', 'Marker', 'none');
    xlabel('Number of Clusters (k)');
    ylabel('Gap Statistic');
    title('Gap Values for Different k');
    grid on;
end
