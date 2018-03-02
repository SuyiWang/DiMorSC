function [G, vert] = swc2graph(filename, clr, flip, trans)
%% swc reader
    if nargin == 0
        filename = 'OP_9.swc';
%         filename = 'test.swc';
%         filename = 'Disperse.swc';
        clr = 'r';
        trans = [0 0 0];
        flip = 0;
    elseif nargin <= 2
        flip = 0;
        trans = [0 0 0];
    elseif nargin <= 3
        trans = [0 0 0];
    end
    
    data = load_v3d_neuron_file(filename);
    L = length(data);
    vert = [data(:,3) data(:,4) data(:,5)];
    % plot3(vert(:,1), vert(:,2), vert(:,3),'.');
    edge = [data(:,1) data(:,7)];
    idx = find(data(:,7)<=0);
    edge(idx,:) = [];
    
    if flip
        vert(:,[1,2]) = vert(:,[2,1]);
    end
    vert = transformvert(vert, trans, [1 1 1]);
    
    % to graph simplify
    G = sparse(edge(:,1), edge(:,2), ones(length(edge), 1), length(vert), length(vert));
    % G = SimpComponent(G, vert);
    
    [I J] = find(G);
    edge = [I J];

    
    hold on;
    patch('faces', edge, 'vertices', vert, 'edgecolor', clr, 'LineWidth', 1);
    % plot3(vert(1,1),vert(1,2),vert(1,3),'r*','markersize',5);
    hold off;
    
    G = G;
    %% translate
%     tt = Tree2SWCtt(G, [vert(:, 1:3) ones(size(vert(:,3)))], 1, 1);
%     save_v3d_swc_file(tt, 'output.swc');