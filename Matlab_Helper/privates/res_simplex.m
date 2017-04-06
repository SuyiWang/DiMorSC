function res_simplex(hFig, dir_path, filename, trans)
% dir_path = '~/Develop/density/Results/Partha/';

if nargin < 4
    trans = [0 0 0];
end

set(hFig, 'Position', [100 100 1200 800])

fp = fopen([dir_path filename],'r');
m = fread(fp, 1, 'int32');
vertex = fread(fp, [4 m], 'double')';
n = fread(fp, 1, 'int32');
edges = fread(fp, [2 n], 'int32')';
k = fread(fp, 1, 'int32');
tris = fread(fp, [3 k], 'int32');
fclose(fp);

vertex(:, [1,2]) = vertex(:, [2,1]);
g = sparse(edges(:, 1) + 1, edges(:, 2) + 1, ones(n, 1), m, m);
DrawGraph(g, vertex(:, 1:3), 'r', 1);