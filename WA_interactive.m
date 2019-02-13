function [] = WA_interactive(V,F,E)

%% Paths
addpath('decimator');
addpath('matlab');

%% Compute Euclidean embedding
if ~exist('E','var')
  [E] = WA_precompute(V,F);
end

%% Pick the anchors
A = pickQuad(V,F);
A = [A repmat([0.33 0.33],size(A,1),1)]; 
B = barycenter(V,F);
fourP = B(A(:,1),:);

%% Generate a grid of weights
[X, Y] = meshgrid(0:0.1:1, 0:0.1:1);
X = reshape(X,size(X,1)*size(X,2),1);
Y = reshape(Y,size(Y,1)*size(Y,2),1);
%%
eye_matrix = eye(4);
W = zeros(size(X,1),4);
for i=1:size(X,1)
  u = X(i);
  v = Y(i);
  W(i,:) = (1-v).*((1-u).*([1 0 0 0]) + u.*([0 1 0 0])) + v.*((1-u).*([0 0 0 1]) + u.*([0 0 1 0]));
end

%% Solve forward
[WA] = WA_forward(E,F,A,W);

%% Convert the barycentric coordinates to point
fid = WA(:,1); 
mask = fid ~= -1;
fid = fid(mask);
bc = [WA(:,2:3) 1.0-sum(WA(:,2:3),2)];
%%
V1 = V(F(fid,1),:);
V2 = V(F(fid,2),:);
V3 = V(F(fid,3),:);
%%
P = repmat(bc(mask,1),1,3) .* V1 + repmat(bc(mask,2),1,3) .* V2 + repmat(bc(mask,3),1,3) .* V3;

%% Plot the points
trisurf(F,V(:,1),V(:,2),V(:,3),'EdgeColor','none'); hold on;
plot3(fourP(:,1),fourP(:,2),fourP(:,3),'+');
plot3(P(:,1),P(:,2),P(:,3),'*'); hold off;
axis equal;
axis vis3d;

% %% Solve inverse for verification
% [W2] = WA_inverse(E,F,A,WA);
% [WA2] = WA_forward(E,F,A,W2);
% 

end