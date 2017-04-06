fprintf('Reading Vertices...\n');
fp = fopen('vert.bin','r');
vert = fread(fp, [4 inf], 'double')';
fclose(fp);


fprintf('Reading Edges...\n');
fp = fopen('edge.bin','r');
edge = fread(fp, [2 inf], 'int32')';
fclose(fp);


fprintf('Processing Edges...\n');
len = length(vert);    
if size(edge,1) ~= 0
    tmpgraph = sparse(edge(:,1) + 1,edge(:,2) + 1,ones(length(edge),1), len, len);
else
    warning('no edges, skipped');
    return;
end
tmpgraph = tmpgraph + tmpgraph';
[I J K] = find(triu(tmpgraph));
edge = [I J];
clear I; clear J; clear K;clear tmpgraph;

fprintf('Reading Triangles...\n');
fp = fopen('triangle.bin','r');
triangles = fread(fp, [3 inf], 'int32')';
fclose(fp);


% fprintf('Reading Tetrahedrons...\n');
% fp = fopen('tetrahedron.bin','r');
% tet = fread(fp, [4 inf], 'int32')';
% fclose(fp);


%%  Write input file for discrete morse
%   edges and triangles are originally using vertex index starting from 1.
fprintf('Writing simplices...\n');
%   since it starts with 0, we put them back to align with original
%   input.
write_output('merge',vert, edge-1, triangles);
% write_tetra([filename 'tet_'], vert, tet - 1);
fprintf('Merged\n');
clear vert; clear edge; clear triangles;

system(['../density/density3D merge.bin inputs/mergevert.txt inputs/mergeedge.txt ' int2str(30)]);
% figure(2);Draw1stable_verbose('inputs/mergevert.txt', 'inputs/mergeedge.txt', [0 0 0]);
Morse_to_tree('inputs/', {'mergevert.txt' 'mergeedge.txt'});