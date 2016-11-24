%%  Plot vertex - edge map
%   input:      Vertex_filename, Edge_filename
%   output:     vert - vertex location and density
%               g - adjacent matrix
%   dependency: DrawGraph.m


function [vert, g] = Draw1stable(vertname, edgename, clr, linesize, flip)
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
    vert = load(vertname);
    edge = load(edgename);
%     sometimes need swap first two dimensions
    if flip
        vert(:,[1,2]) = vert(:,[2,1]);
    end


%%  Create adjacency graph
    g = sparse(edge(:,1),edge(:,2), ones(length(edge),1), length(vert), length(vert));
    DrawGraph(g, vert(:,1:3), clr, linesize);
    
    
%%  Mark critical edges red
    c_edge_idx = find(edge(:,3) == 1);
    edge = edge(c_edge_idx,:);
    
%     g_critical = sparse(edge(:,1),edge(:,2), ones(size(edge,1),1), length(vert), length(vert));
%     DrawGraph(g_critical, vert(:,1:3), 'm',2)
