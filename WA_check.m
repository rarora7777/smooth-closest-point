%% Paths
addpath('decimator');
addpath('matlab');
addpath('toolbox_fast_marching');

%% Read Mesh
[V,F] = readOBJ('data/plane.obj');

%% Compute Euclidean embedding
[E] = WA_precompute(V,F);

%% Some fixed anchors
A = [1 0.3 0.3; 2 0.3 0.3; 3 0.3 0.3];

%% A point on the surface
WA = [4 0.3 0.3];

%% Solve backward
[W] = WA_inverse(E,F,A,WA);

%% Solve forward
[WA2] = WA_forward(E,F,A,W);

%% Check
if max(max(WA - WA2))<0.001
  disp('WA_check OK!');
else
  disl('WA_check FAILED');
end