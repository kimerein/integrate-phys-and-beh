function optimalK = optimalClustersGapStatistic(data, maxK, B)
    % data: NxM data matrix
    % maxK: Maximum number of clusters to consider
    % B: Number of reference datasets for bootstrap sampling

    rng('default');  % For reproducibility

    % Remove rows containing NaN values
    data = data(~any(isnan(data), 2), :);

    % Compute the sum of squared distances for the actual dataset
    actualDists=nan(1,maxK);
    for i=1:maxK
        actualDists(i)=sumofsq(data, i);
    end

    % Generate B reference datasets and compute their distances
    refDists = zeros(B, maxK);
    for i = 1:B
        refData = arrayfun(@(z) min(data(:, z)) + range(data(:, z)) * rand(size(data, 1), 1), 1:size(data, 2), 'UniformOutput', false);
        refData = horzcat(refData{:});
        for k = 1:maxK
            refDists(i, k) = sumofsq(refData, k);
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

function W_allk=sumofsq(data, k)

[idx, C, sumd] = kmeans(data, k, 'EmptyAction', 'singleton', 'Replicates', 10);
uniqueclust=unique(idx);
W_allk=0;
for i=1:length(uniqueclust)
    currclust=uniqueclust(i);
    % sum of pairwise distances for all points in each cluster
    pairwisedist=0;
    for j=1:size(data,1)
        for k=1:size(data,1)
            % if both points in this cluster, add to sum of pairwise
            % distances
            if idx(j)==currclust & idx(k)==currclust
                pairwisedist=pairwisedist+sum((data(j,:)-data(k,:)).^2);
            end
        end
    end
    W_thisk=(1/(2*sqrt(sum(C(currclust,:).^2,'all'))))*pairwisedist;
    W_allk=W_allk+W_thisk;
end

end
