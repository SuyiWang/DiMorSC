function [ T ] = maxspanningtree( G )
    % computes maximum spanning tree for top 5 components in size
    G = -G;
    n = length(G);
%     addpath(genpath('matlab_bgl'));
%     
%     [ci, sizes] = components(G);
%     for i=1:max(ci)
%         sizes(i) = sum(vert(find(ci==i), 4));
%     end
%     [sorted_size, idx] = sort(sizes, 'descend');
%     
%     mark = zeros(max(idx), 1);
%     mark(idx(1)) = ones(1,1);
%     [I J V] = find(G);
%     for i = length(I):-1:1
%         if (mark(ci(I(i)))==0||mark(ci(J(i)))==0)
%             I(i) = []; J(i) = []; V(i) = [];
%         end
%     end
%     G = sparse(I, J, V, n, n)+sparse(J, I, V, n, n);
    G = G + G';
    [I J V] = mst(G);
    T = sparse(I, J, -V, n, n)+sparse(J, I, -V, n, n);
end

