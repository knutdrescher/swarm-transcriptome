function nSwarming = estimateSwarms( mat, swarmfactor, minCells, cellsToConsider)
    
     
     if sum(cellsToConsider)<minCells
         nSwarming = 0;
     else
         if swarmfactor==0
             swarmfactor = 0.01;
         end
         mat_sub = mat(cellsToConsider, :);

         D_angle = pdist( mat_sub(:,3));
         D_angle = min([D_angle; abs(2*pi-D_angle)],[],1);
         D_angle(D_angle>1.4)=100;
         D_orientation = pdist( mat_sub(:,4));
         D_orientation = pi/180*min([D_orientation; abs(180-D_orientation)],[],1);
         D_dist = pdist( mat_sub(:,1:2));

         D =  2*D_dist + D_angle + 0.5*D_orientation;
         Z = linkage(D);

         Y = Z(:,3);

         T = cluster(Z, 'cutoff', 30*swarmfactor*sum(Y)/numel(Y), 'criterion', 'distance');

         newmat = zeros(size(mat_sub,1),5);
         newmat(1:end,1:4) = mat_sub;
         newmat(1:end,5) = T;
         [~,I] = sort(T);
         newmat = newmat(I,:);

         [~, swarmSize] = getClusterMatrix(newmat(:,5), minCells);
         nSwarming = sum(swarmSize);





         
     end



        
end

