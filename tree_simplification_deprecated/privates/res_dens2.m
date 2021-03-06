function res_dens2(hFig, dir_path, filename, zshift)
% dir_path = '~/Develop/density/Results/Partha/';

if nargin < 4
    zshift = 0;
end

set(hFig, 'Position', [100 100 800 600])

fp = fopen([dir_path filename],'r');
m = fread(fp, 1, 'int32');
vertex = fread(fp, [3 m], 'double');
% vertex = fread(fp, [4 m], 'double');
% n = fread(fp, 1, 'int32');
% edges = fread(fp, [2 n], 'int32');
% k = fread(fp, 1, 'int32');
% tris = fread(fp, [3 k], 'int32');
fclose(fp);

% vertex = vertex (:, vertex(3, :)==zshift);

if isempty(vertex)
    disp([filename 'has no vertices, skipped']);
    return;
end

vertex([1,2], :) = vertex([2,1], :);
vertex(1, :) = vertex(1, :) - 1;
vertex(2, :) = vertex(2, :) - 1;
% vertex(3, :) = vertex(3, :);

% cmap = colormap('gray');
% cmap = flipud(cmap);
cmap = colormap('cool');
hold on;
fscatter3(vertex(1,:),vertex(2,:), zshift*ones(size(vertex(1,:))), vertex(3,:)+1, cmap);
axis normal;
cameratoolbar('Show')
% scatter3(vertex(1,:),vertex(2,:), vertex(3,:), 2, vertex(4,:),'filled')

% Draw1stable([dir_path 'outvert8.txt'],[dir_path 'outedge8.txt'],'c', 3);
% 
% edge_info = load('inputs/output_info.txt');
% [~, order] = sort(edge_info(:,1), 'descend');
% for i = 1:length(order)
%     idx = edge_info(order(i), 2) + 1;
%     fac = edges(:, idx) + 1;
%     patch('faces', fac', 'vertices', vertex(1:3, :)', 'edgecolor', 'k', 'linewidth', 5);
%     disp([edge_info(order(i), 1) edge_info(order(i), 3)]);
% end