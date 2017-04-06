% file will write to input/Neocortical_#.bin
% filled triangle

%%  LoadOlfactory_OP(dataset, thd);
dataset = 2; str_data = int2str(dataset);
remove_thd = 32; persist_thd = 40;
keep_top = 30;

% output_config('Results/Neuromuscular_', dataset, remove_thd, persist_thd, keep_top);
% trans = LoadNeuromuscular_blocks(dataset, remove_thd);
% trans = LoadNeuromuscular(dataset, remove_thd);


%%  density3D [input] [output_vertex] [output_edge] [threshold]
% outputname{1} = ['Neuromuscular_' str_data 'vert.txt'];
% outputname{2} = ['Neuromuscular_' str_data 'edge.txt'];
outputname{1} = ['Neuromuscular_'];
outputname{2} = ['Neuromuscular_'];
% system(['../density/density3D inputs/Neuromuscular_2.bin inputs/' outputname{1} ' inputs/' outputname{2} ' ' int2str(persist_thd)]);


%%  Morse_Post(0, [input_folder], [output_filename], [top_#_branches])
% fp = fopen(['supple_data/Neocortical_subset_' int2str(dataset) '_start'],'r');
% if dataset == 1
%     DataNum = 34;
% else
%     DataNum = 21;
% end
% root_info = zeros(DataNum, 3);
% 
% for dataset = 1:DataNum
%     linescan = fgetl(fp);
%     pattern = '(\-?[0-9]+,\-?[0-9]+,\-?[0-9]+)';
%     [start_index, end_index] = regexp(linescan, pattern);
% 
%     root_info(dataset, 1:3) = sscanf(linescan(start_index:end_index),'%d,%d,%d');
% end
% fclose(fp);

%   Must flip to obtain alignment
% root_info(:,[1,2]) = root_info(:,[2,1]);

%%  Draw golden results
% if dataset == 1
%     for branch = 1:DataNum
%     tmpchar = sprintf('%2.2d', branch);
%     swc2graph(['DiademMetric/NC_' sprintf('%2.2d', branch) '.swc'], 'r');
%     end
% else
%     for branch = 1:DataNum
%     tmpchar = sprintf('%2.2d', branch);
%     swc2graph(['DiademMetric/NC_' sprintf('%c', branch - 1 + 'A') '.swc'], 'r');
%     end
% end
% figure(1);
for i = 7:9
    input_filename{1} = [outputname{1} int2str(i) '_vert.txt'];
    input_filename{2} = [outputname{2} int2str(i) '_edge.txt'];
%     text(trans(i,1) + 256, trans(i,2) + 256, trans(i,3), int2str(i),'fontsize', 10);
    figure(2);
    [final_tree, vert, usedidx] = Morse_to_tree('inputs/', input_filename);
%     Draw1stable(['inputs/' input_filename{1}],['inputs/' input_filename{2}], 'c', 2, 0);
end


%%  draw results
figure(1);
% Morse_to_tree('inputs/', {'mergevert.txt' 'mergeedge.txt'});
% clf;
% swc2graph(['DiademMetric/Neocortical_' str_data '.swc'], 'r');
% swc2graph(['inputs/Neocortical_' str_data '.swc'], 'b');
% res_dens(2, 'inputs/',['Neocortical_' str_data '.bin']);
% % swc2graph(['DiademMetric/Olfactory_OP' str_data '.tif_app2.swc'], 'c');
% cameratoolbar('Show')
% colormap cool

