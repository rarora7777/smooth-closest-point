%% Paths
addpath('decimator');
addpath('matlab');

%% Read Mesh
[V,F] = readOFF('data/gargo.off');

%% Compute Euclidean embedding
[E] = WA_precompute(V,F);

%% Launch interactive application
WA_interactive(V,F,E);