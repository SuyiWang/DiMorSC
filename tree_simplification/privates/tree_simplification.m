function tree_simplification(tree, vert, usedidx, output_filename, persist, pos)
    if nargin > 5
        handroot = 1;
    else
        handroot = 0;
    end
    
    if handroot > 0
        [T, handroot] = SimpComponent(tree, vert, pos, usedidx);
    else
        T = SimpComponent(tree, vert);
    end


    %%  Generate files for vaa3D
    disp('writing swc file');
    if handroot > 0
        tt = Tree2SWCtt(T, vert, 2, handroot);
    else
        tt = Tree2SWCtt(T, vert, 2);
    end
    save_v3d_swc_file(tt, ['inputs/' output_filename])
    

    %%  Simplify low persistence branches;
%     if ~is_tree
%         persistence_simplification(persist, 0);
%         Morse_Post(1, 'inputs/', output_filename, persist);
%     end