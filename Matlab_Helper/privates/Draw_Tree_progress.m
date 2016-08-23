function [] = Draw_Tree_progress()
%DRAW_TREE_PROGRESS Summary of this function goes here
%   Detailed explanation goes here
    system(['TreeSimp inputs/tree.bin tree_vert.txt tree_edge.txt ' int2str(1000)]);
    list = load('persistencePairs.txt');
    list = list(list>0);
    for i = length(list):-1:1
        system(['TreeSimp inputs/tree.bin tree_vert.txt tree_edge.txt ' sprintf('%f', list(i) - 1e-6)]);
        Draw1stable('tree_vert.txt','tree_edge.txt','c');
        pause;
    end
end

