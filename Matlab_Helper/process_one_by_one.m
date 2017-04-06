% file will write to input/Olfactory_OP_#.bin
% filled triangle
%%  LoadOlfactory_OP(dataset, thd);
dataset = 8;str_data = int2str(dataset);
remove_thd = 10; persist_thd = 10;
keep_top = 30;
% output_config('Results/OP_', dataset, remove_thd, persist_thd, keep_top);
LoadOlfactory_OP(dataset, remove_thd);
tic

%%  density3D [input] [output_vertex] [output_edge] [threshold]
system(['../density/density3D inputs/Olfactory_OP_' str_data '.bin inputs/outvert.txt inputs/outedge.txt ' num2str(persist_thd)]);
toc
tic;
%%  Morse_Post(0, [input_folder], [output_filename], [top_#_branches])
outname{1} = 'inputs/outvert.txt';
outname{2} = 'inputs/outedge.txt';
Morse_Post(0, outname, ['OP_' str_data '.swc'], keep_top);
toc;
%%  draw results
figure(2);
clf;
swc2graph(['DiademMetric/OP_' str_data '.swc'], 'r');
swc2graph(['inputs/OP_' str_data '.swc'], 'b');
% res_dens(2, 'inputs/',['Olfactory_OP_' str_data '.bin']);
% swc2graph(['DiademMetric/Olfactory_OP' str_data '.tif_app2.swc'], 'c');
cameratoolbar('Show')
colormap cool

%%  compare result
java_cmd = ['java -jar DiademMetric/DiademMetric.jar -G DiademMetric/OP_' str_data '.swc '];
% system([java_cmd '-T DiademMetric/Disperse_' str_data '.swc -D 5 -m false']);
system([java_cmd '-T inputs/OP_' str_data '.swc -D 5 -m false']);
system([java_cmd '-T DiademMetric/Olfactory_OP' str_data '.tif_app2.swc -D 5 -m false']);
system([java_cmd '-T DiademMetric/Olfactory_OP' str_data '.tif_app22.swc -D 5 -m false']);
