function [] = persistence_simplification(depth, progress)
    if nargin == 0
        depth = 30;
        progress = 0;
    end
    system(['./density3D inputs/tree.bin tree_vert.txt tree_edge.txt ' int2str(1000)]);
    list = load('PersistenceValues.txt');
    list = list(list>0);
    
    if progress
        for i = 1:depth
            val = list(length(list) - i);
            system(['./density3D inputs/tree.bin tree_vert.txt tree_edge.txt ' sprintf('%.6f', val - 1e-3)]);
            res_dens();
%             Draw1stable('tree_vert.txt','tree_edge.txt','b',2, false);
            % saveas(1, ['inputs/pic_' int2str(i)], 'fig');
            % saveas(1, ['pic_' int2str(i)], 'png');
        end
    else
        val = list(length(list) - depth);
        system(['./density3D inputs/tree.bin tree_vert.txt tree_edge.txt ' sprintf('%.6f', val - 1e-3)]);
%         Draw1stable('tree_vert.txt','tree_edge.txt','b',2, false);
    end
end

