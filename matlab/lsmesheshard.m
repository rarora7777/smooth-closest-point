function [ NS ] = lsmesheshard( V,F, IDS, S)
% Implementation of Least-squares Meshes, Olga Sorkine and Daniel Cohen-Or
% Input:
% V,F: mesh description
% IDS: #IDSx1 indices of the fixed vertices in S
% S:   #Vxn coordinates of the points to smooth (usually = V)
%
% Output:
% NS:  #Vxn new coordinates 

%display('LSMESHES HARD');

%if (~exist('IDS','var'))
%   IDS = samplemesh(V,F,100); 
%end

%if (~exist('S','var'))
%    S = V;
%end

%% Build left part
L = cotmatrix(V,F);
M = massmatrix(V,F,'voronoi');

%% Build energy matrix
A = L*(M\L);
B = zeros(size(A,1),1);

%% Minimize with hard contraints
NS = min_quad_with_fixed(A,B,IDS,S(IDS,:));

end

