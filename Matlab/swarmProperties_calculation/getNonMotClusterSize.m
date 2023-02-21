function [swarmSize, numberOfSwarms, swarms_colormap_allCells] = getNonMotClusterSize( mat,  minCells, cellsToConsider)
    
     
     if sum(cellsToConsider)<minCells
         swarmSize = [];
         numberOfSwarms = 0;
         swarms_colormap_allCells = [];
     else
     
         mat_sub = mat(cellsToConsider, :);

         D =  pdist( mat_sub);
         Z = linkage(D);

         T = cluster(Z, 'cutoff', 10, 'criterion', 'distance');

         newmat = zeros(size(mat_sub,1),3);
         newmat(1:end,1:2) = mat_sub;
         newmat(1:end,3) = T;
         [~,I] = sort(T);
         newmat = newmat(I,:);

         [swarms_colormap, swarmSize] = getClusterMatrix(newmat(:,3), minCells);
         swarms_colormap(I,:) = swarms_colormap;
         numberOfSwarms = length(swarmSize);

         swarms_colormap_allCells = 0.5*ones(size(mat,1),3);
         swarms_colormap_allCells(cellsToConsider,:) = swarms_colormap;
         swarms_colormap_allCells(2:end+1,:) = swarms_colormap_allCells;
         swarms_colormap_allCells(1,:) = [0, 0, 0];
     
     end
end

