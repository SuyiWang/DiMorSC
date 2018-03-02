%%  Enforce tree structure using shortest path tree
%	Input: Graph file (output from DiMorSC)
%	Output: Tree represented by adjacency list
%	Parameters: (need at least 1 parameters)
%		input_prefix: program will read [input_prefix]_vert.txt and [input_prefix]_edge.txt
%		saddle_threshold: Removes saddles whose density is lower than threshold
%		pos: Root location


function [final_tree, vert] = Morse2sp(pos, input_prefix, saddle_threshold, shift)

		%%  Preset value 
		scale = [1 1 1];
    if nargin < 2
        input_prefix = 'output/'
    end
    if nargin < 3
        saddle_threshold = 0;
    end
    if nargin < 4
        shift = [0 0 0];
    end
		
		%%  Set dependency paths
    %  addpath('vaa3d_matlab_io');
    addpath('privates');
    addpath('matlab_bgl');


    %%  Read graph file
    disp('reading graph file');
    %   figure(2);
    %   dir is also included in the input_prefix
    input_filename = {[input_prefix '_vert.txt'], [input_prefix '_edge.txt']};
    [vert, G] = Draw1stable(input_filename{1}, input_filename{2}, ...
                'c', 0, false, saddle_threshold);
    
    if isempty(G)
        fprintf('%s is empty', input_prefix);
        final_tree = [];
        vert = [];
        return
    end
		%  translate vertex if specified
    vert = transformvert(vert, shift, scale);
    
		%     [I J K] = find(G>0);
		%     K = 2/(vert(I, 4) + vert(J, 4));
		%     K(K<1) = 1;
		%     G=sparse(I, J, K, size(vert, 1), size(vert, 1));
    G(G>0) = 1;
    
		%%  Keep the largest connected component
    disp('remove redundant components');
    newG = G+G';
    [ci, sizes] = components(newG);
    disp('# of components:');
    disp(max(ci));
    
    [sorted_size, idx] = sort(sizes, 'descend');
    % keepmark = find(sorted_size>40);
    keepmark = 1;
    mark = zeros(max(idx), 1);
    mark(idx(keepmark)) = 1;
    [I J V] = find(newG);
    idx = (mark(ci(I)) | mark(ci(J)))==0;
    I(idx) = []; J(idx) = []; V(idx) = [];

    newG = sparse(I, J, V, size(vert, 1), size(vert, 1)); %+sparse(J, I, V, n, n);
    
		
		%%  Find root location
    [G, vert] = compressTree(newG, vert);
    %   DrawGraph(G, vert(:,1:3),'r',1);
    root = findelement(vert, pos);
		
		
		%%  Compute shortest path tree
    [d, pred] = shortest_paths(G, root);
    I = 1:length(d);
    J = pred(I);
    K = d(I);
    I(root) = []; J(root) = []; K(root) = [];
    final_tree = sparse(I, J, K, size(vert, 1), size(vert, 1));
    vert(:, 5) = vert(:, 4);
    vert(:, 4) = d;
    DrawGraph(final_tree, [vert(:, 1:3) d],'c',0.5);

    if (~isTree(final_tree))
        warning('detected: output is not a Tree!');
    end

