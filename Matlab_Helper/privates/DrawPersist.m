function vert = DrawPersist(p_name, b_name, r)

pairs = load(p_name);
% idx1 = find(pairs(:,1) == 1);
% idx2 = find(pairs(:,1) == 2 & pairs(:, 2) > 0);
% idx = [idx1;idx2];
% 
% Epos = pairs(idx,4) + 1;
% Eval = pairs(idx,2);
% 
% fp = fopen(b_name,'r');
% m = fread(fp, 1, 'int32');
% vert = fread(fp, [4 m], 'double')';
% n = fread(fp, 1, 'int32');
% edges = fread(fp, [2 n], 'int32')';
% fclose(fp);
% 
disp(size(pairs, 1))
% 
% MAXC = 64;
% I = edges(Epos, 1);
% J = edges(Epos, 2);
% K = ones(size(Epos, 1), 1);
% K(Eval > 0) = 2;
% K =  Eval;
% maxx = max(K);
% minn = min(K);
% K = fix((K - minn)/(maxx-minn) * (MAXC - 1)) + 1;
% col = colormap(cool(MAXC));
%% DrawPath([25547 25531], vert(:, 1:3), 'c') - 114799
% col = [0 0 0; 1 0 0];
% 
% uK = unique(K);
% L = zeros(length(uK), 1);
% for iuK = 1:length(uK)
%     i = I(K==uK(iuK));
%     j = J(K==uK(iuK));
%     L(iuK) = length(i);
%     patch('faces', [i+1 j+1], 'vertices', vert(:, 1:3), 'edgecolor', col(uK(iuK), :), 'linewidth', r);
% end

% patch('faces', [I+1 J+1], 'vertices', vert(:, 1:3), 'FaceVertexCData', vert(:,4), 'EdgeColor', 'interp', 'LineWidth', r);