function curve = generateCurveOnMesh(V, ER, F, numC, numKP, dist, angle, seed)
    [~, fIdx, bary] = sampleControlPoints(V, F, numC, dist, angle, seed);
    
    numPts = 500;
    [fIdx, bary] = splineOnSurface(V, ER, F, fIdx, bary, 3, numPts+1);
    
    P = getCartesianFromBarycentric(V, F, fIdx, bary);
    
    tsurf(F, V, 'FaceAlpha', 1, 'EdgeAlpha', 0)
    daspect([1 1 1])
    hold on
    line(P(:, 1), P(:, 2), P(:, 3), 'Color', 'r', 'LineWidth', 2)
    
    stride = numPts / numKP;
    curve = struct(...
        'KPI', floor(1:stride:501),...
        'FI', fIdx,...
        'B', bary, ...
        'seed', seed, ...
        'nc', numC, ...
        'dist', dist, ...
        'angle', angle);
    
    scatter3(P(curve.KPI, 1), P(curve.KPI, 2), P(curve.KPI, 3), 50, 'k', 'filled');
    hold off;
    
end
