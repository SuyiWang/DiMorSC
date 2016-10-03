function tt = Tree2SWCtt(Tree, vert, seg_weight, Handroot)

[I J K] = find(Tree);
if nargin == 3
    [maxx, Handroot] = max(vert(:,4));
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
newvert = vert(realpoint,1:4);
newvert(:, 4) = newvert(:, 4)/maxx;

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
% seg_weight only affect the width of segments in vaa3d
tt(:,6) = seg_weight;

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

%% generate distance function
func = zeros(size(newvert,1), 1);
func(newroot) = 0;
for i = 2:length(idx)
    curr = idx(i);
    ancestor = pred(idx(i));
    if (ancestor==-1) 
        continue;
    end
    func(curr) = func(ancestor) + newvert(curr, 4) * sqrt(sum((newvert(curr,:) - newvert(ancestor,:)).*(newvert(curr,:) - newvert(ancestor,:))));
end
% figure;
% scatter3(newvert(:,1), newvert(:,2), newvert(:,3), 3, func, 'filled');
newTree = triu(newTree);
[I J] = find(newTree);
%% Write tree to bin file
    fp = fopen('inputs/tree.bin','w');
    fwrite(fp, length(newvert), 'int32');
    fwrite(fp, [newvert(:, 1:3) func]', 'float');
    fwrite(fp, length(I), 'int32');
    fwrite(fp, [I-1 J-1]', 'int32');
    fwrite(fp, 0, 'int32');
    fclose(fp);
