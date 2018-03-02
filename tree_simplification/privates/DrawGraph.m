function [ output_args ] = DrawGraph( graph, vert, clr, r )
%DRAWGRAPH Summary of this function goes here
%   Detailed explanation goes here
    [gi, gj, gk] = find(graph);    
    if (size(vert, 2) > 3)
        % patch('faces', [gi gj], 'vertices', vert(:, 1:3), 'FaceVertexCData', gk, 'linewidth',r); 
        patch('faces', [gi gj], 'vertices', vert(:, 1:3), 'FaceVertexCData', vert(:,4), 'EdgeColor', 'interp', 'LineWidth', r);
    else
        patch('faces', [gi gj], 'vertices', vert, 'edgecolor', clr, 'linewidth',r); 
    end
    cameratoolbar('Show')
    colorbar;
end

