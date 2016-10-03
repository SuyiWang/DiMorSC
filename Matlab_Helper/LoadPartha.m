%%  This file loads data from Allen Institute
%   Input:          path to folder
%   Output:         tif file, fits file, txt file
%   Requirement:    create a folder named 'inputs' in current directory
%   Dependency:     TriangulationNew.m
%
%   Other output files:    mapinput.txt vert.txt edge.txt triangle.txt


%%  Read file info from input folder
clear all;
addpath('privates');
path = '/media/My Passport/NeuronData/PMD1228/';
% path = '/media/My Passport/NeuronData/PMD1232/';
% path = '/media/My Passport/NeuronData/PMD1463/';
dirpath = dir(path);
fnames = {dirpath.name};


%%  Loop over files
counter = 0;
first = 1;
start = 1; len = 253;
for k = start:(start+len)
    cnt = sprintf('%4.4d', k);
    try
%%      Search for particular type of input file
%         pattern = ['PMD1229&1228-F\d+-\d{4}.\d{2}.\d{2}-\d{2}.\d{2}.\d{2}_PMD1228_\d{1}_' cnt '.jp2'];
        pattern = ['PMD1228_reduce4_' cnt '.tif'];    % 1 - 253
%         pattern = ['PMD1232_reduce4_' cnt '.tif'];      % 1 - 243
%         pattern = ['PMD1463_reduce4_' cnt '.tif'];      % 1 - 299
        fileindex = regexp(fnames,pattern);
        if sum(~cellfun('isempty',fileindex))>1
            warning('More than one records found');
        end
        selectindex = find(~cellfun('isempty',fileindex));
        
        
        %%  Show progress
        lPrompt = 1; % use this for a licensed version
        %lPrompt = 7; %use this for a trial version
    %     dispstr = sprintf('Progress = %d / %d', k, start+len);
        dispstr = fnames{selectindex};
        if (k == start)
            disp(dispstr);
        else
    %   char(8) is the ascii character for "backspace"
    %   dispay the require number of them
            disp([char(8)*ones(1,lStr+lPrompt), dispstr]);        
        end
        lStr = length(dispstr);
    %   End of Show progress
        

        data = imread([path fnames{selectindex}]);
    catch
%%      if a file does not exist, ignore it
        warning('%d does not exist\n', k);
        continue;
    end
    
%%  Flip data, if necessary. This used to align our output with the input.
    data = flip(data, 1);
    
    
%%  All 3-dimension has same data, we pick one of them, truncate if necessary  
    writedata = uint8(data(:,:,1)/256) * 10;
%     image(uint8(data(:,:,1)/256)*10)
    
    
%%  Append data    
    counter = counter + 1;
    if counter==1
        imgdata = zeros([size(writedata) len]);  
        imgdata(:,:, counter) = writedata;
        imwrite(writedata,['inputs/Partha.tif']);
    else
        imgdata(:,:, counter) = writedata;
        imwrite(writedata,['inputs/Partha.tif'],'WriteMode','append');
    end
end
clear data; clear writedata;

PreTriangulation(imgdata);
clear imgdata;

Triangulate('inputs/Partha', 0);
