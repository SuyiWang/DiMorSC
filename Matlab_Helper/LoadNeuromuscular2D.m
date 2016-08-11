%%  This file loads Neuromucsular Layer 1
%   Input:          path - default path
%   Output:         tif file, fits file, txt file.
%   Requirement:    create a folder named 'inputs' in current directory
%   Dependency:     TriangulationNew.m
%
%   Other files:    mapinput.txt vert.txt edge.txt triangle.txt


%%  Data sets - Change this if necessary
%   we have to clean everything otherwise it will run out of memory
    clear all;
    addpath('privates');
%   Neuromucsular data.
    default_path = '/media/My Passport/NeuronData/DIADEM/Neuromuscular Projection Fibers/Subset 2/Image Stacks- 16 bit Gray/';
%     default_path = 'H:/NeuronData/DIADEM/Neuromuscular Projection Fibers/Subset 2/Image Stacks- 16 bit Gray/';
    
%%  Read Translation information
    fp = fopen('supple_data/Neuromuscular_subset_2','r');
    DataNum = 156;
    trans_info = zeros(DataNum, 2);
    for dataset = 1:DataNum
        linescan = fgetl(fp);
        pattern = ['(\-?[0-9]+,\-?[0-9]+)'];
        [start_index, end_index] = regexp(linescan, pattern);
        
        trans_info(dataset, 1:2) = sscanf(linescan(start_index:end_index),'%d,%d');
    end
    fclose(fp);
    global_trans = min(trans_info);
    trans_info(:,1) = trans_info(:,1) - global_trans(1);
    trans_info(:,2) = trans_info(:,2) - global_trans(2);

    
%%  Loop over several data sets (if there is any)
%   The data set contains 6 parts, they should be merged later.
    for dataset = 1:DataNum % Necrotical
        disp(['Processing dataset ' int2str(dataset)]);
        folder_name = sprintf('%3.3d', dataset);
%       Necrotical Layer data.
        path = [default_path folder_name '/'];
        len = dir(path);
        counter = 0;
        thd = 0;
        x = trans_info(dataset, 1);
        y = trans_info(dataset, 2);


%%      Loop over all files of one dataset
        for k = 1:length(len)-2
%%          Show progress
            lPrompt = 1; % use this for a licensed version
            %lPrompt = 7; %use this for a trial version
            dispstr = sprintf('Progress = %d / %d', k, length(len)-2);
            if (k == 1)
                disp(dispstr);
            else
%               char(8) is the ascii character for "backspace"
%               dispay the require number of them
                disp([char(8)*ones(1,lStr+lPrompt), dispstr]);        
            end
            lStr = length(dispstr);
%           End of Show progress
            
            
            cnt = sprintf('%3.3d', k);
            correspond(k) = counter;

            
%%          try data with different format - 1. xxx010.tif; 2. xxx10.tif
            try
                data = imread([path cnt '.TIF']);
            catch
%                 warning('%d does not exist\n', k);
                cnt = sprintf('%2.2d', k);
                correspond(k) = counter;
                try 
                    data = imread([path cnt '.TIF']);
                catch
                    warning('%d REALLY does not exist\n', k);
                    continue;
                end
            end
            
            
%%          Threshold data
%           we shrink data 1) to reduce noise, 2) to save memory
            data = uint8(data/256);
            data(data < 35) = 0;
%             image(uint16(data));
            
%%           On first file, create imgdata to hold the output.
            counter = counter + 1;
            if (dataset == 1 && counter==1)
%               Note 1: subset 1 -> 60; subset 2 -> 85, check the z-value
%               BELOW!
%               Note 2: The data need fliped to be aligned
%               Note 3: in Matlab Vertical is the 1-st dimension. so x-y
%                       needs flipping
                imgdata = zeros([size(data)+max(trans_info(1:DataNum, 1:2))]);
                sumdata = double(data);
            elseif counter == 1
                sumdata = double(data);
            else
                sumdata = sumdata + double(data);
            end
        end
        sumdata = sumdata / 50;
        imgdata(x + 1:x + size(sumdata', 1),...
                y + 1:y + size(sumdata', 2)) = sumdata';
    end
%%  Write tif data
    disp('Writing tif data...');
    imwrite(imgdata(:,:),['inputs/Neuromucsular.tif']);

    
%%  Create simplicial complex
%     disp('Creating triangluation...');
%     
%     TriangulationNew(imgdata, 'inputs/Neuromucsular');
%     disp('***************DONE********************');


