function [pair_fid,pair_bc] = picking(varargin)
  % PICKING demo
  %
  % picking(tsh)
  % picking(V,F)
  %

  pair_fid = [];
  pair_bc  = [];
  
  V = varargin{1};
  F = varargin{2};
  % append phony z-coords
  if size(V,2) ~= 3
    V = [V 0*V(:,1)];
  end
  f = figure();
  tsh = trisurf(F,V(:,1),V(:,2),V(:,3));

  set(tsh,'ButtonDownFcn',@onmeshdown, ...
    'FaceColor','flat', ...
    'FaceLighting','phong', ...
    'CData',0*V(:,1));
  colormap([[255,228,58]/255;[80,64,255]/255]);
  caxis([0 1]);

  % living variables
  down_pos = [];
  middle_down_pos = [];
  hith = [];

  % Update color of mesh depending on whether ray from camera through mouse
  % intersects any faces
  function update()
    % Get current point on front and back clipping planes
    down_pos = get(gca,'currentpoint');
    % Cast ray from front point to back point and intersect with mesh
    [flag,t,lambda] = ...
      ray_mesh_intersect_slow(down_pos(1,:),down_pos(2,:)-down_pos(1,:),V,F);
    % Show all hits in yellow 
    set(tsh,'CData',flag*1.0);
    % mark vertex of closest hit
    flag = flag & t>=0;
    
    candidates        = find(flag);
    [~,id] = min(t(flag));
    
    pair_fid = [pair_fid;candidates(id)];
    pair_bc = [pair_bc;lambda(candidates(id),:)];
    
    if length(pair_fid) == 4
      close(f);
    end
    
  end

  function onmeshdown(src,ev)
    set(gcf,'windowbuttonmotionfcn',@onmeshdrag);
    set(gcf,'windowbuttonupfcn',    @onmeshup);
    update();
  end

  function onmeshdrag(src,ev)
    update();
  end

  function onmeshup(src,ev)
    set(gcf,'windowbuttonmotionfcn',[]);
    set(gcf,'windowbuttonupfcn',    []);
  end
  
  uiwait(f);
end
