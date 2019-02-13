function [W] = WA_inverse(E,F,A,WA)
% [W] = WA_backward(E,F,A,WA);
%
% Input:
% E: #Vx8: 8-dimensional embedding, the metric used is the Euclidean
%    distance in it
% F: #Fx3: mesh faces
% A: #Ax3: anchor points, each row [fid p1 p2], where fid is the face id and p1 and p2 are the first two barycentric coordinates
% WA: #WAx3: point for which we want to find the weights, one row  in the format [fid, p1, p2]
% 
% Output:
% W: #WA x #P: matrix of weights, one row for each requested point

%% Prepare the data to call the WA binary

% Out file
out_path       = [pwd '/t_O.txt'];%[tempname];

% Mesh
mesh_path      = [pwd '/t_M.obj']; %[tempname '.obj'];
writeOBJ(mesh_path,E(:,1:3),F);

% Embedding
embedding_path = [mesh_path '.emb'];
dlmwrite(embedding_path, E, 'delimiter', ' ', 'precision', 12);

% Anchors
anchors_path   = [pwd '/t_A.txt'];%[tempname];
dlmwrite(anchors_path, [A(:,1)-1 A(:,2:end)], 'delimiter', ' ', 'precision', 12);

% Points
points_path   = [pwd '/t_P.txt'];%[tempname];
dlmwrite(points_path, [WA(:,1)-1 WA(:,2:end)], 'delimiter', ' ', 'precision', 12);

%% Call the WA binary
WA_binary = [pwd '/WA'];
command_string = [WA_binary ' i ' mesh_path ' ' anchors_path ' ' points_path ' ' out_path];
unix(command_string);
%command_string

%% Read output
W = dlmread(out_path);

%% Clean up
delete(out_path);
delete(mesh_path);
delete(embedding_path);
delete(anchors_path);
delete(points_path);


end

