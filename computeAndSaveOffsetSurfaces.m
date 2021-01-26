function computeAndSaveOffsetSurfaces(modelName)
    [V, F] = readOBJ("./data/" + modelName + ".obj");
    
    offset = @(V) 0.05*max(1, norm(max(V)-min(V)));
    
    [V1, F1] = signed_distance_isosurface(V, F, 'SignedDistanceType', 'pseudonormal', 'Level', offset(V));
    
    % Just take the largest (#pts) component, presumably that's the one 
    % covering (V,F). Could also find the geometrically largest component.
    [~, CF] = connected_components(F1);
    F1 = F1(CF==1, :);
    [V1, IV1] = remove_unreferenced(V1, F1);
    F1 = IV1(F1);
    
    writeOBJ("./data/" + modelName + "_out_cgal.obj", V1, F1);
    [V2, F2] = signed_distance_isosurface(V, F, 'SignedDistanceType', 'pseudonormal', 'Level', -offset(V));
    
    if numel(F2)==0
        disp('In-mesh is empty.');
    else
        writeOBJ("./data/" + modelName + "_in_cgal.obj", V2, F2);
    end
    
    disp('Saved offset surfaces. Use TetWild to improve mesh quality, and then run computeAndSaveEmbedding');
end
