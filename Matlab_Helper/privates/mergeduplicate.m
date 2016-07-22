function [ G, vert, mergemap, idx] = mergeduplicate( G, vert )
    %% search points with brute force
    nbidx = rangesearch(vert, vert, 1e-3);
    mergemap = 1:size(vert,1);
    keeplist = ones(size(vert,1), 1);
    
    for i = 1 : size(vert, 1)
        if mergemap(i) ~= i
            found = mergemap(i);
        else
            found = 0;
        end
        
        if found>0
            keeplist(i) = 0;
            target1 = find(G(i,:));
            target2 = find(G(:,i));
            G = G + sparse(ones(length(target1), 1)*found, target1, ones(length(target1), 1), size(vert,1), size(vert,1));
            G = G + sparse(target2, ones(length(target2), 1)*found, ones(length(target2), 1), size(vert,1), size(vert,1));
        else
            if size(nbidx{i},2)>1
                mergemap(nbidx{i}) = i;
            end
        end
    end
    %% generate graph
    idx = find(keeplist);
    vert = vert(idx, :);
    G = G(idx, idx);
    save('Graph','G','vert');
end

