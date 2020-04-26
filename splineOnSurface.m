function [path] = splineOnSurface(V, ER, F, intFaces, lambda, order, numSplineEvalPts)

    % (V, F) is the mesh in 3D
    % ER is the n-dimensional embedding of ER
    % intFaces gives the face ids and lambda the barycentric coordinates
    
    % ensure all face idx are valid (all rays intersect the mesh)
    assert(all(intFaces > 0));
    numControlPts = numel(intFaces);
    assert(numControlPts > order);
    
    % Compute a b-spline basis matrix
    knot = [zeros(1, order-1),...
        linspace(0, 1, numControlPts - order + 2),...
        ones(1, order-1)];
    mat = bspline_basismatrix(order, knot, linspace(0, 1, numSplineEvalPts));
    
    WA = WA_forward(ER(1:size(V, 1), :), F, ...
        [intFaces, lambda(:, 1:2)], ...
        mat);
    
    fIdx = WA(:, 1);
    bary = [WA(:, 2) WA(:, 3) 1-WA(:, 2)-WA(:, 3)];


%%
%     % initialize the geodesic computation library
%     global geodesic_library;
%     geodesic_library = 'geodesic_release';
%     mesh = geodesic_new_mesh(V, F);
%     algorithm = geodesic_new_algorithm(mesh, 'exact');
%     
%     path = cell(1, 0);
    
    i = numSplineEvalPts;
    fv = [V(F(fIdx(i), 1), :); V(F(fIdx(i), 2), :); V(F(fIdx(i), 3), :)];
    
    destination = geodesic_create_surface_point('face', fIdx(i), bary(i, :) * fv);
    evalPath = cell(numSplineEvalPts, 1);
    evalPath{i} = destination;
    
    for i = numSplineEvalPts:-1:2
%         fv = [V(F(fIdx(i), 1), :); V(F(fIdx(i), 2), :); V(F(fIdx(i), 3), :)];
%         source_points = {geodesic_create_surface_point('face', fIdx(i), bary(i, :) * fv)};
        
        fv = [V(F(fIdx(i-1), 1), :); V(F(fIdx(i-1), 2), :); V(F(fIdx(i-1), 3), :)];
        destination = geodesic_create_surface_point('face', fIdx(i-1), bary(i-1, :) * fv);
        evalPath{i-1} = destination;
        
%         geodesic_propagate(algorithm, source_points);
%         segment = geodesic_trace_back(algorithm, destination);
%         path = [segment(2:end); path(:)];
    end
%     path = [{destination}; path];


    path = evalPath;
    
    geodesic_delete;
    
    
end
