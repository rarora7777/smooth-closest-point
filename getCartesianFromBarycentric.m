function P = getCartesianFromBarycentric(V, F, fIdx, bary)

    nP = size(fIdx, 1);
    assert(size(bary, 1) == nP);
    
    if size(bary, 2) == size(F, 2)-1
        bary = [bary; 1-sum(bary, 2)];
    end
    
    P = zeros(nP, size(V, 2));
    
    for i=1:nP
        fv = [V(F(fIdx(i), 1), :); V(F(fIdx(i), 2), :); V(F(fIdx(i), 3), :)];
        P(i, :) = bary(i, :) * fv;
    end
end
