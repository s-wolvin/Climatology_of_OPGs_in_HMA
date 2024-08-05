function [clusters,cluster_counts, sum_distance, centroid] = relabel_clusters(clusters, sum_distance, centroid)
% Takes as input a single vector describing the cluster partitioning.
% Outputs a re-labeled clustering in which the clusters are labeled in
% descending order of number of members. Also returns a lookup table for
% the membership counts, original labels, and resulting labels.
 
%% Sort clusters by size to use only the large ones
% create array to hold labels and count
cluster_counts = NaN(max(clusters),2);
cluster_counts(:,1) = 1:max(clusters);

% count number of facets in each original cluster
for j = 1:max(clusters)
   cluster_counts(j,2) = sum(clusters==j); 
end
cluster_counts = flipud(sortrows(cluster_counts,2));

%% Relabel in descending order
% add additional column to list new labels
cluster_counts = [cluster_counts (1:size(cluster_counts,1))'];

% relabel
clusters_orig_labels = clusters;
for j = 1:size(cluster_counts,1)
   clusters(clusters_orig_labels==cluster_counts(j,1)) = cluster_counts(j,3); 
end
cluster_counts = cluster_counts(find(cluster_counts(:,2)>0),:);% keep only nonzero clusters (some initial labelings may screw this up otherwise)

%order sum_distance
sum_distance = sum_distance(cluster_counts(:,1));
centroid = centroid(cluster_counts(:,1),:);

return