function edgeList = s_branch(init_list, vert, verbose)
    edgeList = cell(1, length(init_list));
    for i = 1:length(init_list)
%         disp(i)
        idxlist = init_list{i};
        if isempty(idxlist)
            continue;
        end
        keeplist = zeros(size(idxlist));
        keeplist(1) = 1; keeplist(end) = 1;
        
        curr = 1; next = 2;
        for j = 3:length(idxlist)
            if within_limit(idxlist(curr:j), vert)
                next = j;
            else
                keeplist(next) = 1;
                curr = next; next = j;
            end
        end
        edgeList{i} = idxlist(find(keeplist));
        if verbose
            clf;
            DrawPath(init_list{i}, vert, 'b',3);
            DrawPath(edgeList{i}, vert, 'r',1);
        end
    end