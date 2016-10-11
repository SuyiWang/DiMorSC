%%  This file loads data from Allen Institute
%   Input:          path to folder
%   Output:         tif file, fits file, txt file
%   Requirement:    create a folder named 'inputs' in current directory
%   Dependency:     TriangulationNew.m
%
%   Other output files:    mapinput.txt vert.txt edge.txt triangle.txt


function LoadAllen(path)
%%  Read file info from input folder
clear all;
addpath('privates');
addpath('matlab_bgl');
if nargin == 0
    path = uigetdir('~');
    if path==0
        path = '/media/My Passport/NeuronData/OSU_AllenAAVs/100141454';
        % path = '/media/My Passport/NeuronData/OSU_AllenAAVs/100141563';
        % path = '/media/My Passport/NeuronData/OSU_AllenAAVs/100141780';
    end
end
dirpath = dir(path);
fnames = {dirpath.name};


%%  Loop over files
counter = 0;
first = 1;
start = 191; len = 330 - 191;
%   start = 614; len = 753 - 614;
%   start = 311; len = 449 - 311;
for k = start:(start+len)
    %%  ***ATTENTION*** You might need change the for loop parameter above
    %%  ***ATTENTION*** You might need change the pattern for searching file below.
    cnt = sprintf('%3.3d', k);
    try
        %%  Search for particular type of input file
        pattern = ['102139' cnt '-projection.png'];
        %   pattern = ['102141' cnt '-projection.png'];
        %   pattern = ['102152' cnt '-projection.png'];
        fileindex = regexp(fnames,pattern);
        if sum(~cellfun('isempty',fileindex))>1
            warning('More than one records found');
        end
        selectindex = find(~cellfun('isempty',fileindex));

        data = imread([path '/' fnames{selectindex}]);
    catch
        %%  if a file does not exist, ignore it
        warning('%d does not exist\n', k);
        continue;
    end
    
    
    %%  Show progress
    lPrompt = 1; % use this for a licensed version
    %   lPrompt = 7; % use this for a trial version
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
    
    
    %%  All 3-dimension has same data, we pick one of them, truncate if necessary  
    writedata = data(:, :, 1);
    writedata = rmv_boundary(writedata, [50 70]);

    
    %%  Flip data, if necessary. This used to align our output with the input.
    %   writedata = flip(writedata, 1);
    
        
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


PreTriangulation(imgdata, 10);
clear imgdata;
Triangulate(['inputs/Allen'], 0);
disp('***************DONE********************');