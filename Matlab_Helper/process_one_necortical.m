% file will write to input/Neocortical_#.bin
% filled triangle

%%  LoadOlfactory_OP(dataset, thd);
addpath('privates');
addpath('matlab_bgl');
addpath('vaa3d_matlab_io');
dataset = 1; str_data = int2str(dataset);
remove_thd = 10; persist_thd = 30;
keep_top = 30;

output_config('Results/Neocortical_filter_', dataset, remove_thd, persist_thd, keep_top);
trans = LoadNeocortical_blocks(dataset, remove_thd);


%%  density3D [input] [output_vertex] [output_edge] [threshold]
% outputname{1} = ['Neocortical_filter' str_data 'vert.txt'];
% outputname{2} = ['Neocortical_filter' str_data 'edge.txt'];

outputname{1} = ['Neocortical_'];
outputname{2} = ['Neocortical_'];
% system(['../density/density3D merge.bin inputs/' outputname{1} ' inputs/' outputname{2} ' ' num2str(persist_thd)]);


%%  Draw result
for i = 1:6
    input_filename{1} = [outputname{1} int2str(i) '_vert.txt'];
    input_filename{2} = [outputname{2} int2str(i) '_edge.txt'];
    [vert, g] = Draw1stable(['inputs/' input_filename{1}],['inputs/' input_filename{2}], 'r', 1, 0);
%     clf;
%     maxx = max(vert(:, 4));
%     DrawGraph(g, vert(:,1:3), 'r', 1);
%     Draw1stable_verbose(['inputs/' input_filename{1}],['inputs/' input_filename{2}],[0 0 0]);
%     res_dens(1, '/media/My Passport/NeuronData/Temp/',['NeoC_' int2str(i) '.bin'], [0 0 0]);
    
end
[vert, g] = Draw1stable('inputs/mergevert.txt', 'inputs/mergeedge.txt','r', 2,0);


%%  Morse_Post(0, [input_folder], [output_filename], [top_#_branches])
fp = fopen(['supple_data/Neocortical_subset_' int2str(dataset) '_start'],'r');
if dataset == 1
    % neuron count
    DataNum = 34;
else
    DataNum = 21;
end
root_info = zeros(DataNum, 3);

figure(2);
for branch = 1:DataNum
    linescan = fgetl(fp);
    pattern = '(\-?[0-9]+,\-?[0-9]+,\-?[0-9]+)';
    [start_index, end_index] = regexp(linescan, pattern);

    root_info(branch, 1:3) = sscanf(linescan(start_index:end_index),'%d,%d,%d');
end
fclose(fp);

%   Must flip to obtain alignment
root_info(:,[1,2]) = root_info(:,[2,1]);


%%  Draw standard/trusted results
if dataset == 1
    for branch = 1:DataNum
%         clf;
        DrawGraph(g, vert(:, 1:3), 'c', 1.2);
        swc2graph(['DiademMetric/NC_' sprintf('%2.2d', branch) '.swc'], 'r');
%         pause;
%         saveas(2, ['pics/NeoComp_' int2str(branch)], 'fig');
%         saveas(2, ['pics/NeoComp_' int2str(branch)], 'png');
    end
else
    for branch = 1:DataNum
        tmpchar = sprintf('%2.2d', branch);
        swc2graph(['DiademMetric/NC_' sprintf('%c', branch - 1 + 'A') '.swc'], 'r');
    end
end
% DrawGraph(final_tree, vert(:, 1:3), 'b', 1);
% swc2graph('inputs/Neocortical_1.tif_x33_y203_z31_app2.swc', 'c',1,trans);
% res_dens(1, 'inputs/',['Neocortical_' str_data '.bin'], trans);


%%  draw results
% figure(2);
% clf;
% swc2graph(['DiademMetric/Neocortical_' str_data '.swc'], 'r');
% swc2graph(['inputs/Neocortical_' str_data '.swc'], 'b');
% res_dens(2, 'inputs/',['Neocortical_' str_data '.bin']);
% % swc2graph(['DiademMetric/Olfactory_OP' str_data '.tif_app2.swc'], 'c');
% cameratoolbar('Show')
% colormap cool

%%  compare result
% java_cmd = ['java -jar DiademMetric/DiademMetric.jar -G DiademMetric/OP_' str_data '.swc '];
% system([java_cmd '-T DiademMetric/Disperse_' str_data '.swc -D 5 -m false']);
% system([java_cmd '-T inputs/OP_' str_data '.swc -D 5 -m false']);
% system([java_cmd '-T DiademMetric/Olfactory_OP' str_data '.tif_app2.swc -D 5 -m false']);
% system([java_cmd '-T DiademMetric/Olfactory_OP' str_data '.tif_app22.swc -D 5 -m false']);
