% file will write to input/Neocortical_#.bin
% filled triangle

%%  LoadOlfactory_OP(dataset, thd);
dataset = 1; str_data = int2str(dataset);
remove_thd = 20; persist_thd = 10;
keep_top = 30;

% output_config('Results/Neocortical_', dataset, remove_thd, persist_thd, keep_top);
% trans = LoadNeocortical(dataset, remove_thd);


%%  density3D [input] [output_vertex] [output_edge] [threshold]
% system(['../density/density3D inputs/Neocortical_' str_data '.bin inputs/outvert.txt inputs/outedge.txt ' int2str(persist_thd)]);


%%  Morse_Post(0, [input_folder], [output_filename], [top_#_branches])
fp = fopen(['supple_data/Neocortical_subset_' int2str(dataset) '_start'],'r');
if dataset == 1
    DataNum = 34;
else
    DataNum = 21;
end
root_info = zeros(DataNum, 3);

for dataset = 1:DataNum
    linescan = fgetl(fp);
    pattern = '(\-?[0-9]+,\-?[0-9]+,\-?[0-9]+)';
    [start_index, end_index] = regexp(linescan, pattern);

    root_info(dataset, 1:3) = sscanf(linescan(start_index:end_index),'%d,%d,%d');
end
fclose(fp);

for dataset = 1:DataNum
    pos = root_info(dataset, :);
    Morse_Post(0, 'inputs/', ['Neocortical_' str_data '.swc'], keep_top, pos, trans);
end

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
