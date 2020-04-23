function [WA] = WA_forward(E,F,A,W)
% [WA] = WA_forward(E,F,A,W);
%
% Input:
% E: #Vx8: 8-dimensional embedding, the metric used is the Euclidean
%    distance in it
% F: #Fx3: mesh faces
% A: #Ax3: anchor points, each row [fid p1 p2], where fid is the face id and p1 and p2 are the first two barycentric coordinates
% W: #W x #A: matrix of weights, one row for each requested WA computation
%
% Output:
% WA: #Wx3: weighted averages, one row for each set of provided weights in the format [fid, p1, p2]

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

% Weights
weights_path   = [pwd '/t_W.txt'];%[tempname];
dlmwrite(weights_path, W, 'delimiter', ' ', 'precision', 12);

%% Call the WA binary
if (isunix)
    WA_binary = [pwd '/WA'];
else
    WA_binary = [pwd '/WA_ICC.exe'];
end
command_string = [WA_binary ' f ' mesh_path ' ' anchors_path ' ' weights_path ' ' out_path];

unix(command_string);
%command_string

%% Read output
WA = dlmread(out_path);
WA(:,1) = WA(:,1)+1;

%% Clean up
delete(out_path);
delete(mesh_path);
delete(embedding_path);
delete(anchors_path);
delete(weights_path);

end

