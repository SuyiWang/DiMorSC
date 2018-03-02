function edgeList = resample_edge(init_list)
    edgeList = cell(1, length(init_list));
    for i = 1:length(init_list)
        idxlist = init_list{i};
        if isempty(idxlist)
            continue;
        end
        keeplist = zeros(size(idxlist));
        keeplist(1:5:end) = 1;keeplist(end) = 1;
        
        edgeList{i} = idxlist(find(keeplist));
    end