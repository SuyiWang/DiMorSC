%%  This file loads data from Allen Institute
%   Input:          path to folder
%   Output:         tif file, fits file, txt file
%   Requirement:    create a folder named 'inputs' in current directory
%   Dependency:     TriangulationNew.m
%
%   Other output files:    mapinput.txt vert.txt edge.txt triangle.txt


%%  Read file info from input folder
addpath('privates');
addpath('matlab_bgl');
path = '/media/My Passport/NeuronData/OSU_AllenAAVs/100141454/';
% path = '/media/My Passport/NeuronData/OSU_AllenAAVs/100141563/';
% path = '/media/My Passport/NeuronData/OSU_AllenAAVs/100141780/';
dirpath = dir(path);
fnames = {dirpath.name};


%%  Loop over files
counter = 0;
first = 1;
start = 191; len = 330 - 191;
% start = 614; len = 753 - 614;
% start = 311; len = 449 - 311;
for k = start:(start+len)
%%  Show progress
    lPrompt = 1; % use this for a licensed version
    %lPrompt = 7; %use this for a trial version
    dispstr = sprintf('Progress = %d / %d', k, start+len);
    if (k == start)
        disp(dispstr);
    else
%   char(8) is the ascii character for "backspace"
%   dispay the require number of them
        disp([char(8)*ones(1,lStr+lPrompt), dispstr]);        
    end
    lStr = length(dispstr);
%   End of Show progress


    cnt = sprintf('%3.3d', k);
    try
%%      Search for particular type of input file
        pattern = ['102139' cnt '-projection.png'];
%         pattern = ['102141' cnt '-projection.png'];
%         pattern = ['102152' cnt '-projection.png'];
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
    writedata = data(:, :, 1);
    writedata = rmv_boundary(writedata, [50 70]);
%     subplot(2, 1, 1)
%         image(writedata);
%     subplot(2, 1, 2)
%         image(data(:,:,1));
%%  Flip data, if necessary. This used to align our output with the input.
    writedata = flip(writedata, 1);
    
        
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


PreTriangulation(imgdata);
clear imgdata;
Triangulate(['inputs/Allen' int2str(dataset)], 0);
disp('***************DONE********************');