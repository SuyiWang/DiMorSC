%%  This file loads Climbing Fibers
%   Input:          path - default path
%   Output:         tif file, fits file, txt file.
%   Requirement:    create a folder named 'inputs' in current directory
%   Dependency:     TriangulationNew.m
%
%   Other files:    mapinput.txt vert.txt edge.txt triangle.txt


%%  Data sets - Change this if necessary
    addpath('privates');
%   Climbing Fibers
    default_path = '/media/My Passport/NeuronData/DIADEM/Cerebellar Climbing Fibers/Image Stacks/CF_';


%%  Loop over several data sets (if there is any)
    for dataset = 1:3 % Olfactory data
%       olfactory data
        path = [default_path int2str(dataset) '/'];


        len = dir(path);
        counter = 0;
        thd = 0;

        
%%      Thresholding low-value pixels - Using average intensity
%         for k = 1:length(len)-2
%             disp(counter)
%             cnt = sprintf('%d', k);
%             try
%                 data = imread([path cnt '.tif']);
%             catch
%                 fprintf('%d does not exist\n', k);
%                 cnt = sprintf('%2.2d', k);
%                 try 
%                     data = imread([path cnt '.tif']);
%                 catch
%                     fprintf('%d REALLY does not exist\n', k);
%                     continue;
%                 end
%             end
%             %data = data(139:310,125:456,:);
%             thd = thd + sum(sum(data));
%         end
%         thd = thd/(512*512*length(len)-2);


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
            
            
            cnt = sprintf('%2.2d', k);
            % cnt = sprintf('%3.3d', k);

            
%%          try data with different format - 1. xxx_010.tif; 2. xxx10.tif
            try
                data = imread([path cnt '.tif']);
            catch
                warning('%d does not exist\n', k);
                cnt = sprintf('%d', k);
                try 
                    data = imread([path cnt '.tif']);
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
            
            
%%          Remove low percentage data
%             tmpimage = data;
%             datapercent = double(tmpimage(:,:,1))./double(tmpimage(:,:,1)+tmpimage(:,:,2)+tmpimage(:,:,3));
%             savedata = (datapercent > 0).*double(tmpimage(:,:,1));
%             % savedata = (datapercent < 0.50).*1000;
%             tmpimage(:,:,1) = savedata;        
%             image(uint8(tmp*10));


%%          truncate data, if necessary
            data = sum(data/3, 3);
%             data = flip(data, 1);     % this is used to align output with input.
%%          preview density image
            image(uint8(data/10));
        

%%          write data to 'imgdata' and .tif file.
%           On first file, create imgdata to hold the output.
            counter = counter + 1;
            if counter==1
                imgdata = zeros([size(data) length(len)-2]);  
                imgdata(:,:, counter) = data;
                imwrite(data,['CF_' int2str(dataset) '.tif'])


%%              diamorse input generation
%                 maxcnt = 255;
%                 imwrite(maxcnt-data,'test.pgm');
%                 fp = fopen('test.pgm','a');
%                 fprintf(fp,'\n');
%                 fclose(fp);

            else
                imgdata(:,:, counter) = data;
                imwrite(data,['CF_' int2str(dataset) '.tif'],'WriteMode','append');


%%              diamorse input generation
%                 imwrite(maxcnt-data,'tmp.pgm');
%                 fp = fopen('tmp.pgm','a');
%                 fprintf(fp,'\n');
%                 fclose(fp);
%                 system('cat test.pgm tmp.pgm > test_new.pgm');
%                 system('mv test_new.pgm test.pgm');
            end
        end
        
        
%%      Smooth data (if necessary - esp for DisPerSE)
%         imgdata(:,:,:) = smooth3(imgdata,'gaussian',[3 3 3]);
%         imgdata(imgdata < 10) = 0;
        
        
%%      Write fits file for DisPerSE
%         fitswrite(imgdata,['Olfactory_OP_' int2str(dataset) '.fits']);

        
%%      Create simplicial complex
        clear data; clear writedata;

        PreTriangulation(imgdata);
        clear imgdata;

%         Triangulate(['inputs/Olfactory_OP_' int2str(dataset)], 0);
        Triangulate(['inputs/CF_' int2str(dataset)], 0);
        disp('***************DONE********************');
    end


