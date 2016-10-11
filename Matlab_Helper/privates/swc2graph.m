function output = swc2graph(filename, clr, DrawFlip, HEIGHT)
%% swc reader
    if nargin == 0
        filename = 'OP_9.swc';
%         filename = 'test.swc';
%         filename = 'Disperse.swc';
        DrawFlip = 0;
        clr = 'r';
        HEIGHT = 512;
    elseif nargin <= 2
        DrawFlip = 0;
        HEIGHT = 512;
    end
    
    data = load_v3d_neuron_file(filename);
    L = length(data);
    vert = [data(:,3) data(:,4) data(:,5)];
    % plot3(vert(:,1), vert(:,2), vert(:,3),'.');
    edge = [data(:,1) data(:,7)];
    idx = find(data(:,7)<=0);
    handpick = idx;
    edge(idx,:) = [];
    
    if DrawFlip
        vert(:,2) = HEIGHT - vert(:,2);
    end
    
    % to graph simplify
    G = sparse(edge(:,1), edge(:,2), ones(length(edge), 1), length(vert), length(vert));
    G = SimpComponent(G, vert);
    
    [I J] = find(G);
    edge = [I J];

    
%     hold on;
    patch('faces', edge, 'vertices', vert, 'edgecolor', clr, 'LineWidth', 1);
%     plot3(vert(1,1),vert(1,2),vert(1,3),'r*','markersize',5);
%     hold off;
