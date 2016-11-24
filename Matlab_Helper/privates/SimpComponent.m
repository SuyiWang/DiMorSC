function [ G, handroot ] = SimpComponent( Gin, vertin, pos, realidx )
    handroot = 0;
    if nargin > 2
        existidx = find(realidx > 0);
        mapper(1:length(existidx)) = existidx;
        
        vert(:,1) = vertin(existidx,1) - pos(1);
        vert(:,2) = vertin(existidx,2) - pos(2);
        vert(:,3) = vertin(existidx,3) - pos(3);
        
        dist = sqrt(sum(vert.*vert, 2));
        [sorted_dist, dist_idx] = sort(dist, 'ascend');
        handroot = mapper(dist_idx(1));
    end
    G = Gin + Gin';
    n = length(G);
    % addpath(genpath('matlab_bgl'));
    
    [ci, sizes] = components(G);
    if max(ci) > 1
        disp(max(ci));
        
        if handroot == 0
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
            mark = zeros(max(idx), 1);
            mark(ci(handroot)) = ones(1,1);
            [I J V] = find(G);
            for i = length(I):-1:1
                if (mark(ci(I(i)))==0||mark(ci(J(i)))==0)
                    I(i) = []; J(i) = []; V(i) = [];
                end
            end
            G = sparse(I, J, V, n, n); %+sparse(J, I, V, n, n);
        end
    else
        disp('1');
        G = Gin;
    end
end

