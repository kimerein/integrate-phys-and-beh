function elbowMethod(data, maxK)
    % data: An NxM matrix with N data points and M features.
    % maxK: The maximum number of clusters to evaluate.

    wcss = zeros(1, maxK); % Within-cluster sum of squares

    for k = 1:maxK
        [idx, C, sumd] = kmeans(data, k);
        wcss(k) = sum(sumd);
    end

    % Plotting
    figure;
    plot(1:maxK, wcss, 'o-', 'LineWidth', 2);
    xlabel('Number of Clusters (k)');
    ylabel('Within-Cluster Sum of Squares');
    title('Elbow Method for Optimal k');
    grid on;
end
