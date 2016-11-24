function [G, actualvert] = ToActual(absT, edgeG, edgelist, vert, verbose, abs_idx)
MAXEDGE = sum(cellfun(@length, edgelist));
I = zeros(MAXEDGE, 1);
J = I;
top = 1;

[treeI, treeJ, ~] = find(absT);

for i = 1:length(treeI)
    edgeidx = edgeG(treeI(i), treeJ(i));
    edgelength = length(edgelist{edgeidx});
    if (edgelength<2)
        warning('caught edge length < 2');
    end
    I(top:top+edgelength-2) = edgelist{edgeidx}(1:end-1);
    J(top:top+edgelength-2) = edgelist{edgeidx}(2:end);
    
    if verbose
        DrawPath([treeI(i) treeJ(i)], vert(abs_idx, :), 'ko');
        DrawPath(edgelist{edgeidx}(1:end), vert, 'k');
    end
    
    top = top + edgelength - 1;
end

I = I(1:top-1);
J = J(1:top-1);
actualvert = zeros(length(vert), 1);
actualvert(I) = 1;actualvert(J) = 1;
existed = find(actualvert > 0);
actualvert(existed) = 1:length(existed);

val = ones(length(I), 1);
G = sparse(I, J, val, length(vert), length(vert));