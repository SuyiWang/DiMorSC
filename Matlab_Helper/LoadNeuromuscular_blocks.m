%%  This file loads Neuromucsular Layer 1
%   Input:          path - default path
%   Output:         tif file, fits file, txt file.
%   Requirement:    create a folder named 'inputs' in current directory
%   Dependency:     TriangulationNew.m
%
%   Other files:    mapinput.txt vert.txt edge.txt triangle.txt

function [trans_info slice_info] = LoadNeuromuscular_blocks(folder, selTHD)
%%  Data sets - Change this if necessary
%   we have to clean everything otherwise it will run out of memory
    addpath('privates');
%   Neuromucsular data.
    default_path = '/media/My Passport/NeuronData/DIADEM/Neuromuscular Projection Fibers/Subset 2/Image Stacks- 16 bit Gray/';
%     default_path = 'H:/NeuronData/DIADEM/Neuromuscular Projection Fibers/Subset 2/Image Stacks- 16 bit Gray/';
    
%%  Read Translation information
    % folder should come here:
    % fp = fopen('supple_data/Neuromuscular_subset_2','r');
    DataNum = 156;
    trans_info = zeros(DataNum, 3);
    slice_info = zeros(DataNum, 1);
    %     for dataset = 1:DataNum
    %         linescan = fgetl(fp);
    %         pattern = ['(\-?[0-9]+,\-?[0-9]+)'];
    %         [start_index, end_index] = regexp(linescan, pattern);
    %         
    %         trans_info(dataset, 1:2) = sscanf(linescan(start_index:end_index),'%d,%d');
    %     end
    % fclose(fp);
    trans_info = load('supple_data/Neuromuscular_s2_mat');
    global_trans = min(trans_info);
    trans_info(:,1) = trans_info(:,1) - global_trans(1);
    trans_info(:,2) = trans_info(:,2) - global_trans(2);
    trans_info(:,3) = trans_info(:,3) - global_trans(3);
    trans = [global_trans 0];

    
%%  Loop over several data sets (if there is any)
%   The data set contains 6 parts, they should be merged later.
    % for dataset = 1:DataNum % Necrotical
    for dataset = [1:156] % Necrotical
        disp(['Processing dataset ' int2str(dataset)]);
        folder_name = sprintf('%3.3d', dataset);
%       Necrotical Layer data.
        path = [default_path folder_name '/'];
        len = dir(path);
        counter = 0;
        thd = 0;
        x = trans_info(dataset, 1);
        y = trans_info(dataset, 2);
        z = trans_info(dataset, 3);
        slice_info(dataset) = length(len) - 2;

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

            
%%          try data with different format - 1. xxx010.tif; 2. xxx10.tif
            try
                data = imread([path cnt '.TIF']);
            catch
%                 warning('%d does not exist\n', k);
                cnt = sprintf('%2.2d', k);
                try 
                    data = imread([path cnt '.TIF']);
                catch
                    warning('%d REALLY does not exist\n', k);
                    continue;
                end
            end
            
%%          Downsampling data - this is usually uncessary
%             tmp = zeros(256, 256);
%             for i = 1:512
%                 for j = 1:512
%                     tmp(floor((i+1)/2),floor((j+1)/2)) = tmp(floor((i+1)/2),floor((j+1)/2)) + data(i,j)/4;
%                 end
%             end
%             data = tmp;
            
            
%%          Threshold data
%           we shrink data 1) to reduce noise, 2) to save memory
%             image(data/10);
            data = data/200;
            fit = fspecial('gaussian', [3 3], 0.75);
            data = imfilter(data, fit);
            data(data<32) = 0;
            data = imfilter(data, fit);
            data(data<32) = 0;
            data = imfilter(data, fit);
            data(data<32) = 0;
            data = imfilter(data, fit);
%             image(uint16(data));
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
                imwrite(data', ['/media/My Passport/NeuronData/DIADEM/Neuromuscular Projection Fibers/Subset 2/' int2str(dataset) '.TIFF'], 'TIFF');
            else
                imgdata(:,:,counter) = data';
                imwrite(data', ['/media/My Passport/NeuronData/DIADEM/Neuromuscular Projection Fibers/Subset 2/' int2str(dataset) '.TIFF'], 'TIFF', 'WriteMode','append');
            end
        end
%%  Create simplicial complex

%         addpath('frangi_filter_version2a');
%         options.BlackWhite=false;
%         options.FrangiScaleRange=[1 1];
%         imgdata = FrangiFilter3D(imgdata,options);

        disp('Creating triangluation...');
        PreTriangulation(imgdata, selTHD);
        clear imgdata;

        Triangulate('inputs/', 0, [x y z], dataset);
        outputname{1} = ['Neuromuscular_' int2str(dataset) '_vert.txt'];
        outputname{2} = ['Neuromuscular_' int2str(dataset) '_edge.txt'];
        system(['../density/density3D inputs/' int2str(dataset) '.bin inputs/' outputname{1} ' inputs/' outputname{2} ' 30']);
    end
    disp('***************Triangulation DONE********************');
%     system(['./Merge_Triangulation ' int2str(DataNum)]);


