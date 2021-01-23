function cavityPts = getPointsInsideCavities(V, F)
    [C, CF] = connected_components(F);
    
    nC = max(C);
    
    cavityPts = zeros(nC, 3);
    
    for i=1:max(C)
        FF = F(CF==i, :);
        VV = V(C==i, :);
        
        maxV = max(VV);
        minV = min(VV);
        
        nP = 10;
        
        sDist = 1;
        idx = 1;
        while sDist(idx)>=0
            [x, y, z] = ndgrid(...
                linspace(minV(1), maxV(1), nP), ...
                linspace(minV(2), maxV(2), nP), ...
                linspace(minV(3), maxV(3), nP));
            samplePts = [x(:), y(:), z(:)];
            
            sDist = signed_distance(samplePts, V, FF);
            
            [~, idx] = min(sDist);
            
            nP = nP*2;
        end
        
        cavityPts(i, :) = samplePts(idx, :);
    end
end