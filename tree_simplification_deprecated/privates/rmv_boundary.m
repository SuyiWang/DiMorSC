function [ output_data ] = rmv_boundary( data, range )
%RMV_BOUNDARY Summary of this function goes here
%   Detailed explanation goes here
    nvert = size(data,1) * size(data,2);
    pb = find( data > range(1) & data < range(2));
    lpb = size(pb, 1);
    [I J] = ind2sub(size(data), pb);
    counter = 0;
    edge = zeros(lpb * 8, 2);
    
    for i = -1:1
        for j = -1:1
            if i==0 && j==0
                continue;
            end
            dest = sub2ind(size(data), I+i, J+j);
            edge(counter*lpb+1 : (counter+1)*lpb, :) = [pb dest];
            counter = counter + 1;
        end
    end
    
    G = sparse(edge(:,1), edge(:,2), ones(size(edge, 1), 1), nvert, nvert);
    G = G + G';
    [ci, sizes] = components(G);
    [~,comp_index] = sort(sizes, 'descend');
    output_data = data;
    for i = 1:7
        output_data(ci == comp_index(i)) = 0;
    end
end

