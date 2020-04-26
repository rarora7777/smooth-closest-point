function computeAndSaveEmbedding(modelName, holes)
    [V, F] = readOBJ(['./data/' modelName '.obj']);
    [V1, F1] = readOBJ(['./data/' modelName '_out.obj']);
    
    if isfile(['./data/' modelName '_in.obj'])
        [V2, F2] = readOBJ(['./data/' modelName '_in.obj']);
    else
        disp('In-offset surface does not exist.');
        V2 = zeros(0, 3);
        F2 = zeros(0, 3);
    end
    
    [ER, VS, T] = embedMeshAndSpace(V, F, V1, F1, V2, F2, 1000, 8, holes);
    
    writeMesh(['./data/' modelName '_tet.txt'], VS, ER, T);
    writeMesh(['./data/' modelName '_tri.txt'], V, ER(1:size(V, 1), :), F);
end
