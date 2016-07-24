function [ G, vert ] = SimpComponent( Gin, vertin )
%SIMPCOMPONENT Summary of this function goes here
%   Detailed explanation goes here
    G = Gin + Gin';
    n = length(G);
    % addpath(genpath('matlab_bgl'));
    
    [ci, sizes] = components(G);
    vert = vertin;
    if max(ci) > 1
        disp(max(ci));
        for i=1:max(ci)
            sizes(i) = length(find(ci==i));
        end
        [sorted_size, idx] = sort(sizes, 'descend');

        mark = zeros(max(idx), 1);
        mark(idx(1)) = ones(1,1);
        [I J V] = find(G);
        for i = length(I):-1:1
            if (mark(ci(I(i)))==0||mark(ci(J(i)))==0)
                I(i) = []; J(i) = []; V(i) = [];
            end
        end
        G = sparse(I, J, V, n, n); %+sparse(J, I, V, n, n);
    else
        disp('1');
        G = Gin;
    end
end

