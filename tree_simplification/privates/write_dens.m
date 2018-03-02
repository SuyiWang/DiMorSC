function write_dens(filename, dens, trans,folder_num)
    index = find(dens > 0);
    m = length(index);
    [I J K] = ind2sub(size(dens), index); %% index start from 1
    I = I + double(trans(1) + 1);
    J = J + double(trans(2) + 1);
    K = K + double(trans(3));
    vert(:,:) = [I J K dens(index)];
    
%%  .txt file: for discrete morse - Mainly use this file
    fp = fopen([filename int2str(folder_num) '_dens.bin'],'w');
    fwrite(fp, m, 'int32');
    fwrite(fp, vert', 'double');
    fclose(fp);
