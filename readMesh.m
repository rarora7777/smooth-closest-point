function [V3D, V, F] = readMesh(filename, dim, dimSimplex)
    % make sure it's a triangle or tet mesh
    assert(dimSimplex==3 || dimSimplex==4);
    
    % Note that this function doesn't perform any additional error checking
    % The file is assumed to exist, be available, and conform to 
    % specifications.
    
    f = fopen(filename, 'r');
    
    
    % first, read in |V| and |F|
    nV = str2double(fgetl(f));
    nF = str2double(fgetl(f));
    
    V = zeros(nV, dim);
    V3D = zeros(nV, 3);
%     F = zeros(nF, dimSimplex);
    
    % for each vertex, read in 3D position in a line, followed by nD
    % position in the next line
    for i=1:nV
        V3D(i, :) = str2num(fgetl(f)); %#ok<*ST2NM>
        V(i, :) = str2num(fgetl(f));
    end
    
    if dimSimplex == 3
        formatSpec = '%d %d %d';
    else
        formatSpec = '%d %d %d %d';
    end
    
    % then read in F
    F = fscanf(f, formatSpec, [nF, dimSimplex]);

    % use 1-based indexing
    F = F+1;
    
    fclose(f);
end
