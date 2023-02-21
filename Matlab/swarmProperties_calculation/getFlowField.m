function [ux,uy] = getFlowField( im1, im2)
% ux, uy: flow field (one value per pixel)
% length: speed (one value per pixel)
% averageLength: average speed 
% uxGrid, uyGrid: flow field on grid (one value per grid pixel)
% uSingle, vSingle: single cell flow field
% ccS: cc of original mask
% cc: cc of grid + mask
% nCells: number of cells counted by mask
%
    
    [ux, uy] = HS(im1, im2);
    
    ux = cutHighValues(ux,150);
    uy = cutHighValues(uy,150);
    
    function matrix = cutHighValues(matrix, threshold)
        med = sum(abs(matrix(:)))/numel(matrix);
        [maxi, index] = max(abs(matrix(:)));
        while (maxi > threshold*med)
            matrix(index) = 0;
            med = sum(abs(matrix(:)))/numel(matrix);
            [maxi,index] = max(abs(matrix(:)));
        end
    end

end

