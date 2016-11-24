%%  For every removed edge by the maximum spanning tree, only remove low value parts
%   starting from the lowest-valued edge, all neighbouring edge with value
%   less than 20% (or any given percentage) will be removed and other edges
%   will be add back to the maximum spanning tree.

function G = makeup(restG, edgeG, edgelist, vert, verbose, abs_idx)
percentage = 0.5;
MAXEDGE = sum(cellfun(@length, edgelist));
I = zeros(MAXEDGE, 1);
J = I;
top = 1;
fieldvalue = vert(:,4);

[treeI, treeJ, ~] = find(restG>1e-6);

for i = 1:length(treeI)
    edgeidx = edgeG(treeI(i), treeJ(i));
    edgelength = length(edgelist{edgeidx});
    if (edgelength<2)
        warning('caught edge length < 2');
    end
    
    edgevert = edgelist{edgeidx};
    [minval, minidx] = min(fieldvalue(edgevert));
    
    % going right
    loopvar = minidx;
    while loopvar<=edgelength && fieldvalue(edgevert(loopvar))< (1.0+percentage)*minval
        loopvar = loopvar+1;
    end
    
    if loopvar<edgelength
        I(top:top+edgelength-loopvar-1) = edgelist{edgeidx}(loopvar:end-1);
        J(top:top+edgelength-loopvar-1) = edgelist{edgeidx}(loopvar+1:end);
        top = top + edgelength - loopvar;
        if verbose
        % DrawPath([treeI(i) treeJ(i)], vert(abs_idx, :), 'ko');
        DrawPath(edgelist{edgeidx}(loopvar:end), vert, 'r');
        end 
    end
    
    
    % going left
    loopvar = minidx;
    while loopvar>0 && fieldvalue(edgevert(loopvar)) < (1.0+percentage)*minval
        loopvar = loopvar-1;
    end
    
    if loopvar>1
        I(top:top+loopvar-2) = edgelist{edgeidx}(1:loopvar-1);
        J(top:top+loopvar-2) = edgelist{edgeidx}(2:loopvar);
        top = top + loopvar - 1;
        if verbose
        % DrawPath([treeI(i) treeJ(i)], vert(abs_idx, :), 'ko');
        DrawPath(edgelist{edgeidx}(1:loopvar), vert, 'r');
        end    
    end
end

I = I(1:top-1);
J = J(1:top-1);
val = ones(length(I), 1);
G = sparse(I, J, val, length(vert), length(vert));