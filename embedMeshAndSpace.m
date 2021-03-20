function [ER, VS, T] = embedMeshAndSpace(V, F, V1, F1, V2, F2, n, d, holePos)
% Compute the quasi-geodesic embedding for the mesh V,F
% Euclidean distance in EI approximate the geodesic over V,F
%
% Input:
% V,F: mesh description (must be connected and manifold)
% V1,F1: outer offset surface description (must be genus 0 and manifold)
% V2,F2: inner offset surface description (must be manifold)
% n: number of samples taken over the mesh [default: 1000]
% d: number of dimensions of the embedding [default: 8]
% holePos: #cavity x 3 list of points [default: computed automatically]
% Typically, you'd want the holes to be inside (V2,F2), or inside (V,F) is
% (V2,F2) is empty.
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

    %% tetmesh the inner and outer regions of space
    Vin = [V; V2];
    Fin = [F; F2+size(V, 1)];
    Vout = [V; V1];
    Fout = [F(:, [2 1 3]); F1+size(V, 1)];
    
    % holePos should contain a point inside each cavity
    % We use tetgen to mesh two regions: the outer shell bounded by the
    % mesh and the outer offset surface and the inner shell bounded by the
    % mesh and the inner offset surface. In both cases, we can just
    % position the holes inside the cavities of the inner offset surface,
    % since those are also inside the mesh.
    % Special case: if inner offset surface is empty, we find cavities
    % inside the mesh instead. The second tetgen call does not need any
    % holes in this case.
    % NOTE that the outer offset surface has been assumed to be genus zero.
    if (~exist('holePos','var') || isempty(holePos))
      if numel(F2)==0
          holePos = getPointsInsideCavities(V, F);
      else
          holePos = getPointsInsideCavities(V2, F2);
      end
    end
    
    disp('Computing shell tet meshes...');

    [VS1, T1] = tetgen(Vout, Fout, 'Flags',sprintf('-Yq1.2a%0.17f',16*avgedge(Vout,Fout)^3/(6*sqrt(2))), 'Holes',holePos);
    
    assert(size(boundary_faces(T1), 1)==size(Fout, 1));

    % the inner offset surface can be empty if the offset distance is too high
    if numel(F2)==0
        [VS2, T2] = tetgen(Vin, Fin, 'Flags',sprintf('-Yq1.2a%0.17f',16*avgedge(Vin,Fin)^3/(6*sqrt(2))));
    else
        [VS2, T2] = tetgen(Vin, Fin, 'Flags',sprintf('-Yq1.2a%0.17f',16*avgedge(Vin,Fin)^3/(6*sqrt(2))), 'Holes',holePos);
    end
    
    assert(size(boundary_faces(T2), 1)==size(Fin, 1));

    % Now combine the two
    VS = [VS1; VS2(size(V, 1)+1:end, :)];
    T2temp = T2;
    T2temp(T2 > size(V, 1)) = T2temp(T2 > size(V, 1)) + (size(VS1, 1) - size(V, 1));
    T = [T1; T2temp];
    

    disp('Decimating...');

    %% Simplify
    [~, ~, M] = decimator(V, F, n, 'edgelength');
    D = zeros(length(M), size(V, 1));


    disp('Computing geodesic distances...');
    
    %% Compute distances between the samples
    for i=1:length(M)
        D(i, :) = geodesicdistance(V, F, M(i));
    end
    D = D(:, M);


    disp('Computing n-D embedding...');
    
    %% Embed in n-d, while trying to preserve the distances between the samples
    D = (D + D') ./ 2;
    opts = statset('MaxIter', 2000);
    E  = mdscale(D, d, 'Criterion', 'metricstress', 'Weights', 1./D.^2, 'Options', opts);


    disp('Interpolating with LS-meshes...');
    
    %% Interpolate with LS Meshes
    d = size(E, 2);
    EI = zeros(size(VS, 1), d);
    EI(M, :) = E;

    if (size(M, 1) == size(EI, 1))
        ER = EI;
    else
        ER = lsmesheshard(VS, T, M, EI);
    end
    
    disp('Done!');
end
