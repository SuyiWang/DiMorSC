function [] = persistence_simplification(depth, progress)
    system(['./density3D inputs/tree.bin tree_vert.txt tree_edge.txt ' int2str(1000)]);
    list = load('PersistenceValues.txt');
    list = list(list>0);
    
    if progress
        for i = 1:depth
            val = list(length(list) - i);
            system(['./density3D inputs/tree.bin tree_vert.txt tree_edge.txt ' sprintf('%.6f', val - 1e-3)]);
            Draw1stable('tree_vert.txt','tree_edge.txt','r',1, false);
        end
    else
        val = list(length(list) - depth);
        system(['./density3D inputs/tree.bin tree_vert.txt tree_edge.txt ' sprintf('%.6f', val - 1e-3)]);
        Draw1stable('tree_vert.txt','tree_edge.txt','r',1, false);
    end
end

