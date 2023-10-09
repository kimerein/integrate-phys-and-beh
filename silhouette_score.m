function [s,s_raw] = silhouette_score(all_glm_coef, labels)

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
% Compute the silhouette score for each data point in X
% X: Data matrix (n x m) where n is the number of data points and m is the number of features
% labels: Cluster labels for each data point (n x 1)

n = size(X, 1);
s = zeros(n, 1);
s_raw = zeros(n,1);

for i = 1:n
    % Calculate average distance to points in the same cluster (a(i))
    same_cluster_idx = find(labels == labels(i));
    same_cluster_idx(same_cluster_idx == i) = [];  % remove the point itself
    a_i = mean(sqrt(sum(bsxfun(@minus, X(i,:), X(same_cluster_idx, :)).^2, 2)));

    % Calculate smallest average distance to points in other clusters (b(i))
    b_i = inf;  % initialize with a large value
    unique_labels = unique(labels);
    for lbl = unique_labels'
        if lbl == labels(i), continue; end
        other_cluster_idx = find(labels == lbl);
        distance_to_other_cluster = mean(sqrt(sum(bsxfun(@minus, X(i,:), X(other_cluster_idx, :)).^2, 2)));
        b_i = min(b_i, distance_to_other_cluster);
    end

    % Compute silhouette score for the i-th data point
    s(i) = (b_i - a_i) / max(a_i, b_i);
    s_raw(i) = (b_i - a_i);
end


end