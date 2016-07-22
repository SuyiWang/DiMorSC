function [ output_args ] = DrawGraph( graph, vert, clr, r )
%DRAWGRAPH Summary of this function goes here
%   Detailed explanation goes here
    [gi gj gk] = find(graph);    
    if (size(vert, 1) > 3)
        patch('faces', [gi gj], 'vertices', vert(:, 1:3), 'edgecolor', clr, 'linewidth',r); 
    else
        patch('faces', [gi gj], 'vertices', vert, 'edgecolor', clr, 'linewidth',r); 
    end
end

