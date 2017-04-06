%%  Plot vertex - edge map
%   input:      Vertex_filename, Edge_filename
%   output:     vert - vertex location and density
%               g - adjacent matrix
%   dependency: DrawGraph.m


function [vert, g] = Draw1stable_verbose(vertname, edgename, trans)
    vert = load(vertname);
    edge = load(edgename);
    vert(:,1) = vert(:,1) - trans(1);
    vert(:,2) = vert(:,2) - trans(2);
    vert(:,3) = vert(:,3) - trans(3);

%%  Create adjacency graph
    g = sparse(edge(:,1),edge(:,2), ones(length(edge),1), length(vert), length(vert));
    
    g = g + g';
    [ci, sizes] = components(g);
    disp('# of components:');
    disp(max(ci));
    disp('remove redundant components');
    
    [sorted_size, idx] = sort(sizes, 'descend');
    keepmark = find(sorted_size>0);
    disp([int2str(length(keepmark)) ' left']);
    
    for i = 1:length(keepmark)
        keepidx = find(ci == idx(i));
        gdraw = g(keepidx, keepidx);
        vertdraw = vert(keepidx, 1:3);

        DrawGraph(gdraw, vertdraw(:, 1:3), rand(1,3), 1);
    end