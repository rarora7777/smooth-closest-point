function [ER, VS, T] = embedMeshAndSpace(V, F, V1, F1, V2, F2, n, d)
% Compute the quasi-geodesic embedding for the mesh V,F
% Euclidean distance in EI approximate the geodesic over V,F
%
% Input:
% V,F: mesh description
% V1,F1: outer offset surface
% V2,F2: inner offset surface
% n: number of samples taken over the mesh [default: 500]
% d: number of dimensions of the embedding [default: 8]
%
% Output:
% ER: #Vxn quasi-geodesic embedding
% VS,T: Tet-mesh of the space b/w (V1,F1) and (V2,F2)

    if (~exist('n','var') || isempty(n))
      n = 1000;
    end

    if (~exist('d','var') || isempty(d))
      d = 8;
    end

    %% tetmesh the inner and outer regions of space
    % ASSUMPTIONS: (V, F) has sphere topology, and mean(V) lies inside
    Vin = [V; V2];
    Fin = [F; F2+size(V, 1)];
    Vout = [V; V1];
    Fout = [F(:, [2 1 3]); F1+size(V, 1)];

    [VS1, T1] = tetgen(Vout, Fout, 'Flags',sprintf('-Yq1.2a%0.17f',16*avgedge(Vout,Fout)^3/(6*sqrt(2))), 'Holes',mean(V));

    % the inner offset surface can be empty if the offset distance is too high
    if numel(F2)==0
        [VS2, T2] = tetgen(Vin, Fin, 'Flags',sprintf('-Yq1.2a%0.17f',16*avgedge(Vin,Fin)^3/(6*sqrt(2))));
    else
        [VS2, T2] = tetgen(Vin, Fin, 'Flags',sprintf('-Yq1.2a%0.17f',16*avgedge(Vin,Fin)^3/(6*sqrt(2))), 'Holes',mean(V));
    end

    % Now combine the two
    VS = [VS1; VS2(size(V, 1)+1:end, :)];
    T2temp = T2;
    T2temp(T2 > size(V, 1)) = T2temp(T2 > size(V, 1)) + (size(VS1, 1) - size(V, 1));
    T = [T1; T2temp];

    %% Simplify
    [~, ~, M] = decimator(V, F, n, 'normaldeviation');
    D = zeros(length(M), size(V, 1));

    %% Compute distances between the samples
    for i=1:length(M)
        D(i, :) = geodesicdistance(V, F, M(i));
    end
    D = D(:, M);

    %% Embed in n-d, while trying to preserve the distances between the samples
    D = (D + D') ./ 2;
    E  = mdscale(D, d, 'Criterion', 'metricstress', 'Weights', 1./D.^2);

    %% Interpolate with LS Meshes
    d = size(E, 2);
    EI = zeros(size(VS, 1), d);
    EI(M, :) = E;

    if (size(M, 1) == size(EI, 1))
        ER = EI;
    else
        ER = lsmesheshard(VS, T, M, EI);
    end
end