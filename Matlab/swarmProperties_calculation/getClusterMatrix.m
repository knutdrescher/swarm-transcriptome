function [colormap_clusters, clusterSize] = getClusterMatrix_16x( cluster , minSize)

    [clusters_unique,~, ids] = unique(cluster);
    clusterSize = [];
    colormap_clusters =[];
    
    for j = 1:length(clusters_unique)
        if sum(ids==clusters_unique(j))>=minSize
            clusterSize = [clusterSize, sum(ids==clusters_unique(j))];
            colormap_clusters(find(ids==clusters_unique(j)),1) = rand(1,1);
            colormap_clusters(find(ids==clusters_unique(j)),2) = rand(1,1);
            colormap_clusters(find(ids==clusters_unique(j)),3) = rand(1,1);
        else
            colormap_clusters(find(ids==clusters_unique(j)),1) = 0.5;
            colormap_clusters(find(ids==clusters_unique(j)),2) = 0.5;
            colormap_clusters(find(ids==clusters_unique(j)),3) = 0.5;
        end
    end


end

