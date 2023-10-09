function [confidence] = bootstrap_kmeans_matched(all_glm_coef, k, n_bootstrap)

endminusinds=5;
smooby=8;
indsdelay=21;
coefs_after_outcome=[];
coefs_after_outcome=[coefs_after_outcome smoothMatrix(all_glm_coef(:,[71*3+1+indsdelay:71*3+1+70-endminusinds]),smooby)];
coefs_after_outcome=[coefs_after_outcome smoothMatrix(all_glm_coef(:,[71*4+1+indsdelay:71*4+1+70-endminusinds])+all_glm_coef(:,[71*5+1+indsdelay:71*5+1+70-endminusinds]),smooby)];
coefs_after_outcome=[coefs_after_outcome smoothMatrix(all_glm_coef(:,[71*3+1+indsdelay:71*3+1+70-endminusinds])+all_glm_coef(:,[71*6+1+indsdelay:71*6+1+70-endminusinds]),smooby)];
coefs_after_outcome=[coefs_after_outcome smoothMatrix(all_glm_coef(:,[71*4+1+indsdelay:71*4+1+70-endminusinds])+all_glm_coef(:,[71*5+1+indsdelay:71*5+1+70-endminusinds])+all_glm_coef(:,[71*7+1+indsdelay:71*7+1+70-endminusinds])+all_glm_coef(:,[71*8+1+indsdelay:71*8+1+70-endminusinds]),smooby)];
temp=coefs_after_outcome;
temp=temp./nansum(temp,2);
temp=temp./nanmax(temp,[],2);
% temp(temp<0)=0; % doesn't affect classification
% tempie=temp; tempie(isnan(tempie))=1/size(tempie,2);
% figure(); imagesc(tempie);
% idx_from_glm=kmeans(temp,2,'Replicates',50);

X=temp;

% Perform bootstrap resampling for k-means clustering with matched clusters
% X: Data matrix (n x m) where n is the number of data points and m is the number of features
% k: Number of clusters
% n_bootstrap: Number of bootstrap samples

n = size(X, 1);

% Initial k-means clustering on the entire dataset
[labels, centroids] = kmeans(X, k);

% Matrix to store cluster assignments from each bootstrap sample
all_labels = zeros(n, n_bootstrap);

for i = 1:n_bootstrap
    % Sample with replacement
    idx = randsample(n, n, true);
    sampled_data = X(idx, :);

    % Apply k-means clustering
    [bs_labels, bs_centroids] = kmeans(sampled_data, k);

    % Match the bootstrap clusters to the original clusters based on centroids
    matched_labels = match_clusters(centroids, bs_centroids, bs_labels);

    % Store the matched cluster labels
    all_labels(idx, i) = matched_labels;
end

% For each data point, count the mode (most frequent) cluster assignment
mode_cluster = mode(all_labels, 2);

% Compute confidence as the frequency of the mode cluster assignment divided by number of bootstrap samples
confidence = sum(all_labels == mode_cluster, 2) / n_bootstrap;
end

function matched_labels = match_clusters(orig_centroids, bs_centroids, bs_labels)

% Match clusters from a bootstrap sample to the original clusters
matched_labels = bs_labels;
for i = 1:size(bs_centroids, 1)
    distances = sqrt(sum((bs_centroids(i, :) - orig_centroids).^2, 2));
    [~, best_match] = min(distances);
    matched_labels(bs_labels == i) = best_match;
end

end
