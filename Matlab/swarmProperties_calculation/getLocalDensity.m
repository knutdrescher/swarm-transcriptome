function densities = getLocalDensity( lx, ly, mask )
%calculated local densities. mask will be cut into parts of size lx x ly
%and for each of those the density will be calculated.
    s1 = size(mask,1);
    s2 = size(mask,2);
    dx = floor(s1/lx);
    dy = floor(s2/ly);
    
    densities = zeros(dx*dy, 1);
    
    for j = 1:dx
        for k = 1:dy
            px = 1+lx*(j-1);
            py = 1+ly*(k-1);
            sub = mask(px:px+lx-1, py:py+ly-1);
            densities(dy*(j-1)+k) = mean(sub(:)); 
        end
    end
            


end

