function [ G ] = Loadswc( filename, pos )
% load swc, set point closest to pos as root, keep only that component
    [G, vert] = swc2graph(filename, 'c');
    [G, handroot] = SimpComponent( G, vert, pos, 1:size(vert,1) );
    G = G + G';
    tt = Tree2SWCtt(G, [vert ones(length(vert), 1)], 2, handroot);
    save_v3d_swc_file(tt, [filename '_reroot.swc']);
end

