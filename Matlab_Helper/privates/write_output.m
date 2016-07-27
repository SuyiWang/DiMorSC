%%  Generate simplicial complex files of different format
%   Input:          filename -- for output
%                   vertex - 4*n matrix, first 3 dimension describes xyz
%                       location, 4-th dimension describes density.
%                   Edge and triangles represented by vertex indices
%
%   Output:         .off file: for data validation
%                   .txt file: for discrete morse
%                   .NDnet: for DisPerSE
%
function write_output(filename, vertex, edges, triangles)


    m = size(vertex, 1);
    n = size(edges, 1);
    k = size(triangles, 1);
%%  .off file: for data validation
%     fp = fopen([filename '.off'],'w');
%     fprintf(fp,'OFF\n');
%     fprintf(fp, '%d %d %d\n', m, k, n);
%     fprintf(fp, '%f %f %f\n', vertex(:,1:3)');
%     fprintf(fp, '3 %d %d %d\n', triangles');
%     fprintf(fp, '2 %d %d\n', edges');
%     fclose(fp);

    
%%  .txt file: for discrete morse - Mainly use this file
    fp = fopen([filename '.bin'],'w');
%     fprintf(fp, '%d\n', m);
%     fprintf(fp, '%f %f %f %f\n', vertex');
%     fprintf(fp, '%d\n', n);
%     fprintf(fp, '%d %d\n', edges');
%     fprintf(fp, '%d\n', k);
%     fprintf(fp, '%d %d %d\n', triangles');
    fwrite(fp, m, 'int32');
    fwrite(fp, vertex', 'double');
    fwrite(fp, n, 'int32');
    fwrite(fp, edges', 'int32');
    fwrite(fp, k, 'int32');
    fwrite(fp, triangles', 'int32');
    fclose(fp);

    
%%  .NDnet: for DisPerSE
%     fp = fopen([filename '.a.NDnet'],'w');
%     fprintf(fp, 'ANDNET\n3\n');
%     fprintf(fp, '%d\n', m);
%     fprintf(fp, '%f %f %f\n', vertex(:,1:3)');
%     fprintf(fp, '2 %d\n', k);
%     fprintf(fp, '%d %d %d\n', triangles');
%     fprintf(fp, '[ADDITIONAL_DATA]\nfield_value\n');
%     fprintf(fp, '0\n');
%     fprintf(fp, '%f\n', vertex(:,4)');
%     fclose(fp);