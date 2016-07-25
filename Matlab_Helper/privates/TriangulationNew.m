%%  Generate triangulation from a given matrix
%   Input:          matrix -- m-by-n-by-k double matrix representing density map
%                   filename -- for output
%
%   Output:         vertex, edge, triangle matrix
%                   Edge and triangles represented by vertex indices
%
%   Requirement:    create a folder named 'input' in current directory
%
%   Dependency:     1) write_output.m
%                   2) ./test -- cpp file compiled from triangulation.cpp
%
%   Other files:    mapinput.txt vert.txt edge.txt triangle.txt


function [vert, edge, triangles] = TriangulationNew(density_map, filename)


%%  Truncate data using a threshold
    fprintf('Smoothing data...\n');
    THD = 20;
    thd_map = smooth3(density_map(:,:,:),'gaussian',[5 5 5], 1); 
    index = find(thd_map > THD);

    real_density_map = smooth3(density_map(:,:,:),'gaussian',[5 5 5]);
    density_map = zeros(size(real_density_map));
    density_map(index) = max(real_density_map(index), 1e-6);

    
%%  Write sparse vertex location file for triangulation
    fprintf('Preparing input...\n');
    len = length(index);
    vert = zeros(len, 4);
    fp = fopen('mapinput.txt','w');
%     for idx = 1:len
%         [i j k] = ind2sub(size(density_map), index(idx));
%         vert(idx,:) = [i j k density_map(i,j,k)];
%     end
    [I J K] = ind2sub(size(density_map), index);
    vert(:,:) = [I J K density_map(index)];
    mnk = size(density_map);
    fprintf(fp, '%d %d %d %d\n', mnk(1), mnk(2), mnk(3), len);
    fprintf(fp, '%d %d %d %f\n', vert');
    fclose(fp);

    
%%  Triangulate - 1) do not fill the cube; 2) fill the cube interier with a tetrahedron
    system('./test');
    % system('./test_fill');

    
%%  Load output
    fprintf('Reading Vertices...\n');
    vert = load('vert.txt');
    fprintf('Reading Edges...\n');
    edge = load('edge.txt');
    
    fprintf('Processing Edges...\n');
    len = length(vert);    
    tmpgraph = sparse(edge(:,1),edge(:,2),ones(length(edge),1), len, len);
    tmpgraph = tmpgraph + tmpgraph';
    [I J K] = find(triu(tmpgraph));
    edge = [I J];
    
    fprintf('Reading Triangles...\n');
    triangles = load('triangle.txt');
    
    
%%  Write input file for discrete morse
%   edges and triangles are originally using vertex index starting from 1.
    fprintf('Writing simplices...\n');
    write_output(filename,vert, edge-1, triangles-1);