function tt = Tree2SWCtt(Tree, vert, Handroot)

[I J K] = find(Tree);
if nargin == 2
    Handroot = 1;
end
n = length(vert);

idx = zeros(n,1);
idx(I) = 1; idx(J) = 1;

realpoint = find(idx);
newroot = find(realpoint == Handroot);

if isempty(newroot)
    newroot = 1;
end

newTree = Tree(realpoint,realpoint);
newvert = vert(realpoint,1:3);

newTree = newTree + newTree';
[~, dt ,~, pred] = dfs(newTree, newroot);

idx = find(pred == 0);
pred(idx) = -1;

[~,idx] = sort(dt);

n = length(newvert);
tt = ones(n,7)* (-1);
tt(:,1) = 1:n;
tt(:,2) = 2;
tt(:, 3:5) = newvert(idx,1:3);
tt(:,6) = 1;

newpred = pred(idx);
reverseidx(idx) = 1:n;
for i=1:length(idx)
    % dt(idx)
    if (newpred(i) ==-1)
        newpred(i) = -1;
    else
       newpred(i) = reverseidx(newpred(i));
    end
end


% pred = pred(idx);
% nonroot = find(pred ~= -1);
% dt(idx) = idx;
% pred(nonroot) = idx(pred(nonroot));

tt(:,7) = newpred;