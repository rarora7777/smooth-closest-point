function writeMesh(filename, V3D, V, F)
    f = fopen(filename, 'w');
    nV = size(V, 1);
    dim = size(V, 2);
    nF = size(F, 1);
    dimSimplex = size(F, 2);
    
    % use 0-based indexing
    F = F-1;
    
    % make sure it's a triangle or tet mesh
    assert(dimSimplex==3 || dimSimplex==4);
    % ensure V3D is actually 3D
    assert(size(V3D, 2) == 3);
    % ensure |V| = |V3D|
    assert(size(V3D, 1) == nV)
    
    % first, write out |V| and |F|
    fprintf(f, "%d\n%d\n", nV, nF);
    
    % for each vertex, write out 3D position ina line, followed by nD
    % position in the next line
    for i=1:nV
        fprintf(f, "%f %f %f\n", V3D(i, 1), V3D(i, 2), V3D(i, 3));
        for j=1:dim-1
            fprintf(f, "%f ", V(i, j));
        end
        fprintf(f, "%f\n", V(i, dim));
    end
    
    % then write out F, one face/cell per line
    for i=1:nF
        if dimSimplex==3
            fprintf(f, "%d %d %d\n", F(i, 1), F(i, 2), F(i, 3));
        elseif dimSimplex==4
            fprintf(f, "%d %d %d %d\n", F(i, 1), F(i, 2), F(i, 3), F(i, 4));
        end
    end
    
    fclose(f);
end
