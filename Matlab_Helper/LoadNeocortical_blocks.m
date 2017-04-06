%%  This file loads Necrotical Layer 1
%   Input:          path - default path
%   Output:         tif file, fits file, txt file.
%   Requirement:    create a folder named 'inputs' in current directory
%   Dependency:     TriangulationNew.m
%
%   Other files:    mapinput.txt vert.txt edge.txt triangle.txt

function trans_info = LoadNeocortical_blocks(folder, selTHD, shift)
    if folder == 1
        DataNum = 6;
    elseif folder == 2
        DataNum = 10;
    end
    
    
    %%  Data sets - Change this if necessary
    addpath('privates');
    %   Necrotical Layer data.
    default_path = ['/media/My Passport/NeuronData/DIADEM/Neocortical Layer 1 Axons/Subset ' int2str(folder) '/Image Stacks/0'];

    
    %%  Read Translation information
    fp = fopen(['supple_data/Neocortical_subset_' int2str(folder)],'r');
    trans_info = zeros(DataNum, 3);
    for dataset = 1:DataNum
        linescan = fgetl(fp);
        pattern = ['(\-?[0-9]+,\-?[0-9]+,\-?[0-9]+)'];
        [start_index, end_index] = regexp(linescan, pattern);
        
        trans_info(dataset, 1:3) = sscanf(linescan(start_index:end_index),'%d,%d,%d');
    end
    fclose(fp);
    global_trans = min(trans_info);
%     trans_info(:,1) = trans_info(:,1) - global_trans(1);
%     trans_info(:,2) = trans_info(:,2) - global_trans(2);
%     trans_info(:,3) = trans_info(:,3) - global_trans(3);
    

    
    %%  Loop over several data sets (if there is any)
    %   The data set contains 6 parts, they should be merged later.
    for dataset = 1:DataNum % Necrotical


        %   Necrotical Layer data.
        path = [default_path int2str(dataset) '/'];
        len = dir(path);
        counter = 0;
        x = trans_info(dataset, 1);
        y = trans_info(dataset, 2);
        z = trans_info(dataset, 3);


        %%  Loop over all files of one dataset
        for k = 1:length(len)-2
            %%  Show progress
            lPrompt = 1; % use this for a licensed version
            %lPrompt = 7; %use this for a trial version
            dispstr = sprintf('Progress = %d / %d', k, length(len)-2);
            if (k == 1)
                disp(dispstr);
            else
            %   char(8) is the ascii character for "backspace"
            %   dispay the require number of them
                disp([char(8)*ones(1,lStr+lPrompt), dispstr]);        
            end
            lStr = length(dispstr);
            %   End of Show progress
            
            
            cnt = sprintf('%d', k);

            
            %%  try data with different format - 1. xxx010.tif; 2. xxx10.tif
            try
                data = imread([path cnt '.tif']);
            catch
                warning('%d does not exist\n', k);
                cnt = sprintf('%2.2d', k);
                try 
                    data = imread([path cnt '.tif']);
                catch
                    warning('%d REALLY does not exist\n', k);
                    continue;
                end
            end
            
            
            %%  Threshold data
            % data(data < 35) = 0;
            % image(uint8(data));
            % data = data(60:240, 31:190);
            fit = fspecial('gaussian', [5 5], 0.75);
            data = imfilter(data, fit);
            data(data<10) = 0;
            data = imfilter(data, fit);
%             data(data<15) = 0;
%             data = imfilter(data, fit);
            
%%          On first file, create imgdata to hold the output.
            counter = counter + 1;
            if (counter==1)
%               Note 1: subset 1 -> 60; subset 2 -> 85, check the z-value
%               BELOW!
%               Note 2: The data need fliped to be aligned
%               Note 3: in Matlab Vertical is the 1-st dimension. so x-y
%                       needs flipping
                imgdata = zeros([size(data') length(len) - 2]);
                imgdata(:,:,counter) = data';
            else
                imgdata(:,:,counter) = data';
            end
        end
        %%  Create simplicial complex
        disp('Creating triangluation...');
        PreTriangulation(imgdata, selTHD);
        clear imgdata;

        Triangulate('/media/My Passport/NeuronData/Temp/NeoC_', 0, [x y z], dataset);
        outputname{1} = ['Neocortical_' int2str(dataset) '_vert.txt'];
        outputname{2} = ['Neocortical_' int2str(dataset) '_edge.txt'];
        system(['../density/density3D /media/My\ Passport/NeuronData/Temp/NeoC_' int2str(dataset) '.bin inputs/' outputname{1} ' inputs/' outputname{2} ' 20']);
    end
    disp('*************** ALL DONE ********************');
