%%  This file loads data from Allen Institute
%   Input:          path to folder
%   Output:         tif file, fits file, txt file
%   Requirement:    create a folder named 'inputs' in current directory
%   Dependency:     TriangulationNew.m
%
%   Other output files:    mapinput.txt vert.txt edge.txt triangle.txt


%%  Read file info from input folder
addpath('privates');
path = '/media/My Passport/NeuronData/OSU_AllenAAVs/100141454/';
dirpath = dir(path);
fnames = {dirpath.name};


%%  Loop over files
counter = 0;
first = 1;
start = 235; len = 10;
for k = start:(start+len)
% for k = 191:330


    disp(counter)
    cnt = sprintf('%3.3d', k);
    try
%%      Search for particular type of input file
%         pattern = ['PMD1229&1228-F\d+-\d{4}.\d{2}.\d{2}-\d{2}.\d{2}.\d{2}_PMD1228_\d{1}_' cnt '.jp2'];
        pattern = ['102139' cnt '-projection.png'];
        fileindex = regexp(fnames,pattern);
        if sum(~cellfun('isempty',fileindex))>1
            warning('More than one records found');
        end
        selectindex = find(~cellfun('isempty',fileindex));

        data = imread([path fnames{selectindex}]);
    catch
%%      if a file does not exist, ignore it
        warning('%d does not exist\n', k);
        continue;
    end
    
    
%%  All 3-dimension has same data, we pick one of them, truncate if necessary  
    writedata = data(200:500,550:800,1);
    
    
%%  Append data    
    counter = counter + 1;
    if counter==1
        imgdata = zeros([size(writedata) len]);  
        imgdata(:,:, counter) = writedata;
        imwrite(writedata,['inputs/allen.tif']);
    else
        imgdata(:,:, counter) = writedata;
        imwrite(writedata,['inputs/allen.tif'],'WriteMode','append');
    end
end

TriangulationNew(imgdata, 'inputs/Allen');
