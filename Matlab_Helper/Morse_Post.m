addpath('vaa3d_matlab_io');
addpath('privates');
addpath('matlab_bgl');

figure;
% [vert, G] = Draw1stable('inputs/outvert.txt','inputs/outedge.txt');
[a b] = VTKread('inputs/Neocortical.fits_c100.up.NDskl.S000.a.vtk');
G = b{1};
vert = [b{2} b{10}];
DrawGraph(G, vert, 'c', 2);

%% adjust position and scale if necessary
trans = [0 0];
scale = [1 1];

% transform to original scale - depend on whether the original data is
% downsampled.
vert = transformvert(vert, trans, scale);
% critical index == 2 -> saddle

% clean with branch
[absG, edgeG, edgeList, abs_idx] = ToAbstract_branch(G, 1, vert);
absT = maxspanningtree(absG);

hold on;
plot3(vert(abs_idx,1),vert(abs_idx,2),vert(abs_idx,3),'r.');
DrawGraph(absT, vert(abs_idx,:),'k',2);
hold off;

G = ToActual(absT, edgeG, edgeList, vert, 0, abs_idx);
G = SimpComponent(G, vert);

restG = logical(absG>1e-6)  - logical(absT>1e-6);
Gadd = makeup(restG, edgeG, edgeList, vert, 0, abs_idx);
% idx = find(b{9}==2 | b{9} ==3);
% hold on;
% plot3(vert(idx,1),vert(idx,2),vert(idx,3),'r.');
% hold off;
% G = saddleclean(G, vert, idx, 1);

DrawGraph(G+Gadd, vert, 'b', 3);
DrawGraph(Gadd, vert, 'k', 2);

tt = Tree2SWCtt(G+Gadd, vert);
save_v3d_swc_file(tt, 'OP_1_disperse.swc')