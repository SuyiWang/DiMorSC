function z_info = image_sticher()
%%  This file finds z-alignment of images

%%  Data sets - Change this if necessary
    addpath('privates');
    default_path = '/media/My Passport/NeuronData/DIADEM/Neuromuscular Projection Fibers/Subset 2/Image Stacks- 16 bit Gray/';
    start_data = 9;
    
%%  Read Translation information
    % folder should come here:
    fp = fopen('supple_data/Neuromuscular_subset_2','r');
    DataNum = 156; % 156 - 2;
    trans_info = zeros(DataNum, 2);
    z_info = zeros(DataNum, 3);
    slice_info = zeros(DataNum, 1);
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
    trans = [global_trans 0];

    h = figure(1);
    clf;
    marked = zeros(DataNum, 1);
%     marked = ones(DataNum, 1) * -1;
%     marked(63) = 0; marked(64) = 0;
%%  Loop over several data sets (if there is any)
    while(sum(marked) < DataNum)
    for dataset = 1:DataNum % Necrotical
        if marked(dataset) == 1 || marked(dataset) == -1
            continue;
        end
        disp(['Processing dataset ' int2str(dataset)]);
        folder_name = sprintf('%3.3d', dataset);
%       Necrotical Layer data.
        path = [default_path folder_name '/'];
        len = dir(path);
        counter = 0;
        thd = 0;
        x = trans_info(dataset, 1);
        y = trans_info(dataset, 2);
        slice_info(dataset) = length(len) - 2;
        
        hasoverlap = 0;
        for dataloop = 1:DataNum
            if (marked(dataloop) == 1) &&...
               (abs(trans_info(dataloop, 1) - x) <= 512 &&...
               abs(trans_info(dataloop, 2) - y) <= 512)
                    bx = min(x, trans_info(dataloop, 1)) + 512;
                    by = min(y, trans_info(dataloop, 2)) + 512;
                    sx = max(x, trans_info(dataloop, 1));
                    sy = max(y, trans_info(dataloop, 2));
                    if (bx - sx)*(by - sy) > (512 * 512 / 20)
                        if trans_info(dataloop, 1) < sx
                            xrng = [sx - trans_info(dataloop, 1) bx - trans_info(dataloop, 1)];
                        else
                            xrng = [1 bx - sx];
                        end
                        
                        if trans_info(dataloop, 2) < sy
                            yrng = [sy - trans_info(dataloop, 2) by - trans_info(dataloop, 2)];
                        else
                            yrng = [1 by - sy];
                        end
                        valid = count_pixel(alldata{dataloop}, xrng, yrng);
                        if valid > 2000
                            hasoverlap = dataloop;
                            break;
                        end
                    end
            end
        end

        if (dataset ~= start_data && hasoverlap <= 0)
            disp(['no overlap skipping' int2str(dataset)]);
            continue;
        end
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
                cnt = sprintf('%2.2d', k);
                try 
                    data = imread([path cnt '.TIF']);
                catch
                    warning('%d REALLY does not exist\n', k);
                    continue;
                end
            end
            
%%          Threshold data
%           we shrink data 1) to reduce noise, 2) to save memory
            data = data/8000;
%             image(uint16(data));
%%          On first file, create imgdata to hold the output.
            counter = counter + 1;
            
            if (counter==1)
%               Note 1: subset 1 -> 60; subset 2 -> 85, check the z-value
%               BELOW!
%               Note 2: The data need fliped to be aligned
%               Note 3: in Matlab Vertical is the 1-st dimension. so x-y
%                       needs flipping
                newdata = false([size(data) length(len)-2]);
                newdata(:,:,counter) = logical(data');
            else
                newdata(:,:,counter) = logical(data');
            end
        end
        
        if dataset == start_data
            alldata = cell(DataNum, 1);
            alldata{dataset} = newdata;
            z_info(dataset, 1) = 0;
            z_info(dataset, 2) = -1;
            z_info(dataset, 3) = 100;
        else
            pdata = alldata{hasoverlap};
            alldata{dataset} = newdata;
            z_info(dataset, 1) = alignimage(pdata, z_info(hasoverlap, 1), newdata, trans_info(hasoverlap, :), trans_info(dataset, :), dataset);
            z_info(dataset, 2) = hasoverlap;
            z_info(dataset, 3) = valid;
        end
        
        disp(['Dataset ' int2str(dataset) ' ' int2str(z_info(dataset))]);
        marked(dataset) = 1;
%%      Draw pic for location processed
        figure(h);
        rectangle('Position',[x, y, 512, 512],'FaceColor',rand(1,3));
    end
    end
    disp('***************All DONE********************');
end

function rtn = alignimage(alldata, zorigin, newdata, trp, trnew, dataset)
    swapped = 1;
    if size(newdata, 3) > size(alldata, 3)
        tmp = newdata;
        newdata = alldata;
        alldata = tmp;
        swapped = -1;
        
        tmp = trp;
        trp = trnew;
        trnew = tmp;
    end
    
    % we are sure that size of all data is longer than newdata.
    minx = min(trp(1), trnew(1)); miny = min(trp(2), trnew(2));
    maxx = max(trp(1), trnew(1)) + 512; maxy = max(trp(2), trnew(2)) + 512;
    
    maxz = size(alldata, 3) - 1;
    minz = size(newdata, 3) - 1;
    
    maxmatch = 0;
    rtn = -1;
    for z = -minz:maxz
%%      Show progress
        lPrompt = 1; % use this for a licensed version
        %lPrompt = 7; %use this for a trial version
        dispstr = sprintf('Progress = %d / %d', z, maxz);
        if (z == -minz)
            disp(dispstr);
        else
%           char(8) is the ascii character for "backspace"
%           dispay the require number of them
            disp([char(8)*ones(1,lStr+lPrompt), dispstr]);        
        end
        lStr = length(dispstr);
%       End of Show progress

        clear expall; clear expnew;
        [expall, expnew] = extendpicture(z, maxx, minx, maxy, miny, maxz,...
                                          minz, trp, trnew, alldata, newdata);
        
        if z == -minz
            topr = uint8(sum(expall,3));
            topg = uint8(sum(expnew,3));
            topb = zeros(size(topr), 'uint8');
            figure(2);
            image(cat(3, topr, topg, topb));
            saveas(2, ['./pics/' int2str(dataset) '_over.png']);
        end
        
        matched = expall & expnew;
        counter = sum(matched(:));
        
        if counter > maxmatch
            maxmatch = counter;
            rtn = z;
%             drawmatch(expall, expnew);
        end
        

    end
    
    [expall, expnew] = extendpicture(rtn, maxx, minx, maxy, miny, maxz,...
                                     minz, trp, trnew, alldata, newdata);
    drawmatch(expall, expnew, dataset);
    
    rtn = zorigin + rtn * swapped;
end

function drawmatch(alldata, newdata, dataset)
    topr = uint8(permute(sum(alldata,2), [1 3 2]));
    topg = uint8(permute(sum(newdata,2), [1 3 2]));
    topb = zeros(size(topr), 'uint8');
    figure(3);
    image(cat(3, topr, topg, topb));
    saveas(3, ['./pics/' int2str(dataset) '_side_y.png']);
    
    topr = uint8(permute(sum(alldata,1), [3 2 1]));
    topg = uint8(permute(sum(newdata,1), [3 2 1]));
    topb = zeros(size(topr), 'uint8');
    image(cat(3, topr, topg, topb));
    saveas(3, ['./pics/' int2str(dataset) '_side_x.png']);
end

function [expall, expnew] = extendpicture(z, maxx, minx, maxy, miny, maxz,...
                                          minz, trp, trnew, alldata, newdata)
    if z < 0
        expall = false([maxx - minx, maxy - miny, maxz - z]);
        expnew = false([maxx - minx, maxy - miny, maxz - z]);

        allx = trp(1) - minx; ally = trp(2) - miny;
        expall((allx + 1):(allx + 512), (ally + 1):(ally + 512), ...
                (-z):(maxz - z)) = alldata;

        newx = trnew(1) - minx; newy = trnew(2) - miny;
        expnew((newx + 1):(newx + 512), (newy + 1):(newy + 512), ...
                1: (minz + 1)) = newdata;
    elseif z + minz > maxz
        expall = false([maxx - minx, maxy - miny, z + minz+1]);
        expnew = false([maxx - minx, maxy - miny, z + minz+1]);

        allx = trp(1) - minx; ally = trp(2) - miny;
        expall((allx + 1):(allx + 512), (ally + 1):(ally + 512), ...
                1:(maxz + 1)) = alldata;

        newx = trnew(1) - minx; newy = trnew(2) - miny;
        expnew((newx + 1):(newx + 512), (newy + 1):(newy + 512), ...
                (z + 1): (minz + z + 1)) = newdata;
    else
        expall = false([maxx - minx, maxy - miny, maxz+1]);
        expnew = false([maxx - minx, maxy - miny, maxz+1]);

        allx = trp(1) - minx; ally = trp(2) - miny;
        expall((allx + 1):(allx + 512), (ally + 1):(ally + 512), ...
                1:(maxz + 1)) = alldata;

        newx = trnew(1) - minx; newy = trnew(2) - miny;
        expnew((newx + 1):(newx + 512), (newy + 1):(newy + 512), ...
                (z + 1): (minz + z + 1)) = newdata;
    end
end

function validpixel = count_pixel(data, xrng, yrng)
    cover_data = data(xrng(1):xrng(2), yrng(1):yrng(2), :);
    validpixel = sum(cover_data(:));
end