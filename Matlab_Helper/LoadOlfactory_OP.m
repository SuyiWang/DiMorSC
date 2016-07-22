%%  This file loads Olfactory, Necrotical Layer 1 and Neuromuscular data
%   Input:          path - default path
%   Output:         tif file, fits file, txt file.
%   Requirement:    create a folder named 'input' in current directory
%   Dependency:     TriangulationNew.m
%
%   Other files:    mapinput.txt vert.txt edge.txt triangle.txt
%   Some notes:
%                   3 5 7 8 div2 does not work - Olfactory


%%  Data sets - Change this if necessary
    addpath('privates');
%   olfactory data
    default_path = '/media/My Passport/NeuronData/DIADEM/Olfactory Projection Fibers/Image Stacks/OP_';
%   Necrotical Layer data.
%   default_path = '/media/My Passport/NeuronData/DIADEM/Neocortical Layer 1 Axons/Subset 1/Image Stacks/0';
%   Neuromuscular data
%   default_path = '/media/My Passport/NeuronData/DIADEM/Neuromuscular Projection Fibers/Image Stacks/0';


%%  Loop over several data sets (if there is any)
    for dataset = 1:1 % Olfactory data
%   for dataset = 1:6 % Necrotical
%   for dataset = 1:22 % Neuromuscular


%       olfactory data
        path = [default_path int2str(dataset) '/'];
%       Necrotical Layer data.
%       path = [default_path int2str(dataset) '/'];
%       Neuromuscular data
%       path_cnt = sprintf('%2.2d', dataset);
%       path = [default_path path_cnt '/'];


        correspond = zeros(500,1);
        len = dir(path);
        counter = 0;
        thd = 0;

        
%%      Thresholding low-value pixels - Using average intensity
%         for k = 1:length(len)-2
%             disp(counter)
%             cnt = sprintf('%d', k);
%             correspond(k) = counter;
%             try
%                 data = imread([path cnt '.tif']);
%             catch
%                 fprintf('%d does not exist\n', k);
%                 cnt = sprintf('%2.2d', k);
%                 correspond(k) = counter;
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
            disp(counter)
            cnt = sprintf('%d', k);
            % cnt = sprintf('%3.3d', k);
            correspond(k) = counter;

            
%%          try data with different format - 1. xxx_010.tif; 2. xxx10.tif
            try
                data = imread([path cnt '.tif']);
            catch
                fprintf('%d does not exist\n', k);
                cnt = sprintf('%2.2d', k);
                correspond(k) = counter;
                try 
                    data = imread([path cnt '.tif']);
                catch
                    fprintf('%d REALLY does not exist\n', k);
                    continue;
                end
            end

            
%%          Remove low percentage data
%             tmpimage = data;
%             datapercent = double(tmpimage(:,:,1))./double(tmpimage(:,:,1)+tmpimage(:,:,2)+tmpimage(:,:,3));
%             savedata = (datapercent > 0).*double(tmpimage(:,:,1));
%             % savedata = (datapercent < 0.50).*1000;
%             tmpimage(:,:,1) = savedata;        
%             image(uint8(tmp*10));


%%          truncate data, if necessary
%             data = data(250:300, 140:180);
%%          preview density image
%             image(uint8(data));
        

%%          write data to 'imgdata' and .tif file.
%           On first file, create imgdata to hold the output.
            counter = counter + 1;
            if counter==1
                imgdata = zeros([size(data) length(len)-2]);  
                imgdata(:,:, counter) = data;
                imwrite(data,['Olfactory_OP' int2str(dataset) '.tif'])
%                 imwrite(data,['Neocortical_' int2str(dataset) '.tif']);
%                 imwrite(data,['Neuromuscular_' int2str(dataset) '.tif'])


%%              diamorse input generation
%                 maxcnt = 255;
%                 imwrite(maxcnt-data,'test.pgm');
%                 fp = fopen('test.pgm','a');
%                 fprintf(fp,'\n');
%                 fclose(fp);

            else
                imgdata(:,:, counter) = data;
                imwrite(data,['Olfactory_OP' int2str(dataset) '.tif'],'WriteMode','append');
%                 imwrite(data,['Neocortical_' int2str(dataset)'.tif'],'WriteMode','append');
%                 imwrite(data,['Neuromuscular_' int2str(dataset) '.tif'],'WriteMode','append');


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
%         imgdata(imgdata<thd*100) = 0;
%%      Write fits file for DisPerSE
%         fitswrite(imgdata,['Neocortical_G3_' int2str(dataset) '.fits']);
%         fitswrite(imgdata,['Neuromuscular_G3_' int2str(dataset) '.fits']);
%%      Create simplicial complex
        TriangulationNew(imgdata, 'inputs/3D');
end


