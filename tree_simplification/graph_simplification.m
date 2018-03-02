function graph_simplification(filename, thd)
    %clf
    load(filename);
    % finalrg / dens
    vert = [pts(:, 1:3) dens(:)];
    % parpool(3);
    [absG, edgeG, edgeList, abs_idx] = ToArc(finalrg, 0, vert, 2);
    % delete(gcp('nocreate'))
    % vert = [vert degG];
    edgeList = s_branch(edgeList, vert, thd, false);
    
    [newG, usedidx] = ToActual(absG, edgeG, edgeList, vert, 0, abs_idx);
    DrawGraph(newG, vert(:, 1:3), 'b', 2);
    print('-dpng', '-r200', filename);
end