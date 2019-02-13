function [ NV, NF, M ] = decimator(V,F,target,type)
  % DECIMATOR 
  %
  % [ NV, NF, M ] = normals(V,F,target)
  %
  % Decimate the mesh using quadric error metric to the number of vertices
  % specified by target
  %
  % Inputs:
  %  V        #V x 3  matrix of vertex coordinates
  %  F        #F x 3  matrix of indices of triangle corners
  %  target   1x1 number of vertices in the target mesh
  %  type     ['quadric': quadric edge collapse priority,
  %            'edgelenght': edge lenght priority,
  %            'normaldeviation': normal deviation priority]
  %
  % Output:
  %  NV       target x 3 matrix of vertex coordinates
  %  NF       #NF x 3    matrix of indices of triangle corners
  %  M        target x 1 maps indices of vertices of NV into V

 if ~exist('target','var')
     target = size(V,1)/2;
 end

 if ~exist('type','var')
     type = 'quadric';
 end

 
 %% Write a temporary file with the mesh
 writeOFF('temp.off',V,F);
 
 %% Decimate the mesh
 if strcmp(type,'quadric')
    system(['decimator\decimator.exe -m temp.off -s ' num2str(target)]);
 elseif strcmp(type,'edgelength')
    system(['decimator\decimator.exe -m temp.off -l -s ' num2str(target)]);
 else
    system(['decimator\decimator.exe -m temp.off -n -s ' num2str(target)]);
 end
 
 %% Read the decimated mesh
 [NV,NF] = readOFF('temp.off-decimated.off');
 
 %% Read the map between decimated vertices and original mesh
 fid = fopen('temp.off-decimated.map', 'r');
 M = fscanf(fid, '%d');
 M = M+1; % indices in matlab starts from 1
 fclose(fid);
 
 %% Clean up the temporary files
 delete('temp.off');
 delete('temp.off-decimated.off');
 delete('temp.off-decimated.map');
end

