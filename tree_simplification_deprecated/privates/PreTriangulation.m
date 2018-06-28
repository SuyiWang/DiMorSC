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


function vert = PreTriangulation(density_map, selTHD, folder_num)
%%  Truncate data using a threshold
    fprintf('Smoothing data...\n');
    % thd_map = density_map;
    thd_map = smooth3(density_map(:,:,:),'gaussian',[13 13 3], 0.98); % used for OP set
%     thd_map = smooth3(density_map(:,:,:),'gaussian',[11 11 11], 0.98);
%       thd_map = smooth3(density_map(:,:,:),'gaussian',[3 3 3], 0.70); %for NeuroMuscular
%       thd_map = density_map; % do not smooth - Neocortical
%       thd_map(thd_map < selTHD) = 0;
%       thd_map = smooth3(thd_map(:,:,:),'gaussian',[3 3 3], 0.70);
        clrTHD = -1e-6;
        thd_map(thd_map < clrTHD) = 1e-6;
    index = find(thd_map > selTHD);
    density_map = thd_map;
    clear thd_map;
    

%     % This blocked is used for first selecting a bigger area, then uses narrower filter.    
%     real_density_map = smooth3(density_map(:,:,:),'gaussian',[7 7 5]);
%     density_map = zeros(size(real_density_map));
%     density_map(index) = max(real_density_map(index), 1e-6);
%     clear real_density_map;
    

%%  Write sparse vertex location file for triangulation
    fprintf('Preparing input...\n');
    len = length(index);
	
	%%  skip data with too many points
%	if len > 3e6
%		fp = fopen('abandoned.txt','a');
%		fprintf(fp, '%d %d\n', folder_num, len);
%		fclose(fp);
%		len = 0
%		index = [];
%	end
	
    vert = zeros(len, 4);
%     fp = fopen('mapinput.txt','w');
    fp = fopen([int2str(folder_num) '_dens.bin'],'w');

    [I J K] = ind2sub(size(density_map), index); %% index start from 1
    mnk = size(density_map);
    vert(:,:) = [I J K density_map(index)];
    clear density_map;
    
    fwrite(fp, [mnk(1) mnk(2) mnk(3) len]', 'int32');
    fwrite(fp, vert', 'double');
%     fprintf(fp, '%d %d %d %d\n', mnk(1), mnk(2), mnk(3), len);
%     fprintf(fp, '%d %d %d %f\n', vert');
    fclose(fp);

