function [ER, VS, T] = embedMeshAndSpace(V, F, V1, F1, V2, F2, n, d, holePos)
% Compute the quasi-geodesic embedding for the mesh V,F
% Euclidean distance in EI approximate the geodesic over V,F
%
% Input:
% V,F: mesh description
% V1,F1: outer offset surface
% V2,F2: inner offset surface
% n: number of samples taken over the mesh [default: 1000]
% d: number of dimensions of the embedding [default: 8]
% holePos: #cavity x 3 list of points [default: mean(V)]
% Typically, you'd want the holes to be inside (V,F) and outside (V2, F2)
%
% Output:
% ER: #Vxn quasi-geodesic embedding
% VS,T: Tet-mesh of the space b/w (V1,F1) and (V2,F2)
%       If (V2,F2) is empty, then (VS,T): tet-mesh b/w (V1,F1) and (V,F)

    if (~exist('n','var') || isempty(n))
      n = 1000;
    end

    if (~exist('d','var') || isempty(d))
      d = 8;
    end
    
    if (~exist('holePos','var') || isempty(holePos))
      holePos = mean(V);
    end

    %% tetmesh the inner and outer regions of space
    Vin = [V; V2];
    Fin = [F; F2+size(V, 1)];
    Vout = [V; V1];
    Fout = [F(:, [2 1 3]); F1+size(V, 1)];

    [VS1, T1] = tetgen(Vout, Fout, 'Flags',sprintf('-Yq1.2a%0.17f',16*avgedge(Vout,Fout)^3/(6*sqrt(2))), 'Holes',holePos);

    % the inner offset surface can be empty if the offset distance is too high
    if numel(F2)==0
        [VS2, T2] = tetgen(Vin, Fin, 'Flags',sprintf('-Yq1.2a%0.17f',16*avgedge(Vin,Fin)^3/(6*sqrt(2))));
    else
        [VS2, T2] = tetgen(Vin, Fin, 'Flags',sprintf('-Yq1.2a%0.17f',16*avgedge(Vin,Fin)^3/(6*sqrt(2))), 'Holes',holePos);
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
    opts = statset('MaxIter', 2000);
    E  = mdscale(D, d, 'Criterion', 'metricstress', 'Weights', 1./D.^2, 'Options', opts);

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
