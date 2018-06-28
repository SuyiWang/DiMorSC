%%  Plot vertex - edge map
%   input:      Vertex_filename, Edge_filename
%   output:     vert - vertex location and density
%               g - adjacent matrix
%   dependency: DrawGraph.m


function [vert, g] = Draw1stable2(vertname, edgename, clr, linesize, flip, zshift)
%%  Set up parameters based on input
    if (nargin < 3)
        % filename = 'test.fits_c50.up.NDskl.a.segs';
        clr = 'c';
        linesize = 2;
        flip = true;
    end
    if (nargin < 4)
        linesize = 2;
        flip = true;
    end
    if (nargin < 5)
        flip = true;
    end
    if (nargin < 6)
        zshift = -1;
    end
%     
%     fp = fopen(filename,'r');
%     fwrite = fopen('tmp.segs','w');
%     
%     sread = fgets(fp);
%     while ~isempty(sread)
%         element = sscanf(sread,'%f');
%         if length(element) > 2
%             fprintf(fwrite,'%s',sread);
%         end
%         
%         sread = fgets(fp);
%     end
    if ~exist(vertname, 'file')
        disp([vertname ' does not exist']);
        vert = [];
        g = [];
        return;
    end
    
    if ~exist(edgename, 'file')
        disp([edgename ' does not exist']);
        vert = [];
        g = [];
        return;
    end
    vert = load(vertname);
    edge = load(edgename);
    
    
    if isempty(edge)
        disp('empty edge, skipping');
        g = [];
        return
    end
    
	
    
%     sometimes need swap first two dimensions
    if flip
        vert(:,[1,2]) = vert(:,[2,1]);
    end
	
	if zshift > -1
		tmp = [vert(:, 1) vert(:, 2) zshift*ones(size(vert, 1), 1) vert(:, 3)];
		vert = tmp;
	end

    %%  Create adjacency graph
    g = sparse(edge(:,1),edge(:,2), ones(size(edge,1),1), length(vert), length(vert));
    if linesize > 0
        DrawGraph(g, vert(:,1:3), clr, linesize);
    end
    
    
%%  Mark critical edges red
%     c_edge_idx = find(edge(:,3) == 1);
%     edge = edge(c_edge_idx,:);
%     
%     g_critical = sparse(edge(:,1),edge(:,2), ones(size(edge,1),1), length(vert), length(vert));
%     DrawGraph(g_critical, vert(:,1:3), 'r',5)
%     c_vert_idx = find(vert(:, 5) == 0);
%     hold on;
%     plot3(vert(c_vert_idx, 1), vert(c_vert_idx, 2), vert(c_vert_idx, 3), 'ro', 'markersize',10);
