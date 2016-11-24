function res_dens(hFig, dir_path, filename)
% dir_path = '~/Develop/density/Results/Partha/';
set(hFig, 'Position', [100 100 1200 800])

fp = fopen([dir_path filename],'r');
m = fread(fp, 1, 'int32');
vertex = fread(fp, [4 m], 'double');
n = fread(fp, 1, 'int32');
edges = fread(fp, [2 n], 'int32');
k = fread(fp, 1, 'int32');
tris = fread(fp, [3 k], 'int32');
fclose(fp);

rand_idx = floor(rand(length(vertex), 1) * 1);
vertex = vertex (:, rand_idx==0);

idx = find(vertex(4,:) > 30);
vertex([1,2], :) = vertex([2,1], :);
vertex = vertex (:, idx);

cmap = colormap('gray');
cmap = flipud(cmap);
hold on;
fscatter3(vertex(1,:),vertex(2,:), vertex(3,:), log(vertex(4,:)+1), cmap);
axis normal;
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