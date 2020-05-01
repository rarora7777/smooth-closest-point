function [vIdx, fIdx, bary] = sampleControlPoints(V, F, numSample, targetDist, targetTheta, randomSeed)

    if (nargin < 6)
        warning('No random seed provided. Results cannot be replicated!');
        rng('shuffle');
    else
        rng(randomSeed, 'twister')
    end
    
    s = rng;
    fprintf('Using random seed %d\n', s.Seed);
    
    nV = size(V, 1);
%     nF = size(F, 1);
    
    % normalize mesh size
    % targetDist is given w.r.t a size-normalized mesh
    bboxSize = max(max(V) - min(V));
    V = V/bboxSize;
    
    % targetTheta is given in degrees. Convert to radian.
    targetTheta = pi/180 * targetTheta;
    
    % compute face centers
%     FC = (V(F(:, 1), :) + V(F(:, 2), :) + V(F(:, 3), :))/3;
    
    % compute unit normal vectors
    N = per_vertex_normals(V, F);
    N = N./vecnorm(N, 2, 2);

    
    vIdx = zeros(numSample, 1);
    fIdx = zeros(numSample, 1);
    bary = zeros(numSample, 3);
    
    % first point is just sampled uniformly
    vIdx(1) = randi(nV);
    
    % initialize library
    global geodesic_library;
    geodesic_library = 'geodesic_release';
    mesh = geodesic_new_mesh(V, F);
    algorithm = geodesic_new_algorithm(mesh, 'exact');

    
    [f, b] = findPositionInFaces(F, vIdx(1));
    fIdx(1) = f;
    bary(1, :) = b;
%     distGeodesicOld = inf(nV, 1);
    
    avoid = false(nV, 1);
    
    for i=2:numSample
        srcIdx = vIdx(i-1);
        source_points = {geodesic_create_surface_point('vertex', srcIdx, V(srcIdx, :))};
        geodesic_propagate(algorithm, source_points);
        
        [~, distGeodesic] = geodesic_distance_and_source(algorithm);
        
        distGauss = acos(dot(N, repmat(N(srcIdx, :), nV, 1), 2));
        
        err = (distGeodesic - targetDist).^2 + ((distGauss - targetTheta)/pi).^2;
        err(avoid) = inf;
        
        [~, minIdx] = min(err);
        vIdx(i) = minIdx;
        [f, b] = findPositionInFaces(F, minIdx);
        fIdx(i) = f;
        bary(i, :) = b;
        avoid(distGeodesic < 0.5*targetDist) = true;
%         distGeodesicOld = distGeodesic;
    end
    

    geodesic_delete;
end

function [f, b] = findPositionInFaces(F, v)
    b = [0 0 0];
    for i=1:3
        f = find(F(:, i) == v, 1);
        if (numel(f)==1)
            b(i) = 1;
            return;
        end
    end
end
