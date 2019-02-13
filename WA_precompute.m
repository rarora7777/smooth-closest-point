function [ER] = WA_precompute(V,F,n,d)
% Compute the quasi-geodesic embedding for the mesh V,F
% Euclidean distance in EI approximate the geodesic over V,F
%
% Input:
% V,F: mesh description
% n: number of samples taken over the mesh [default: 500]
% d: number of dimensions of the embedding [default: 8]
%
% Output:
% ER: #Vxn quasi-geodesic embedding

if (~exist('n','var') || isempty(n))
  n = 500;
end

if (~exist('d','var') || isempty(d))
  d = 8;
end

%% Simplify
[NV,NF,M] = decimator(V,F,n,'normaldeviation');
D = zeros(length(M),size(V,1));

%% Compute distances between the samples
for i=1:length(M)
    D(i,:) = geodesicdistance(V,F,M(i));
end
D = D(:,M);

%% Embed in n-d, while trying to preserve the distances between the samples
T = (D + D') ./ 2;
E  = mdscale(T,d,'Criterion','metricstress','Weights',1./T.^2);

%% Interpolate with LS Meshes
d = size(E,2);
EI = zeros(size(V,1),d);
EI(M,:) = E;

if (size(M,1) == size(EI,1))
    ER = EI;
else
    ER = lsmesheshard(V, F, M, EI);

end

