function [skip] = Triangulate( filename, fill, trans, id )
    skip = 0;
    if (nargin <= 2)
        trans = [0 0 0];
        id = '';       
    else
        % translate the vertices, then write to specific id.
    end
    filename = [filename int2str(id)];
    
    %%  Triangulate - 0) do not fill the cube; 1) fill the cube interier with a tetrahedron
%     system('rm vert.bin');
%     system('rm edge.bin');
%     system('rm triangle.bin');
    fprintf('Trianguating...');
    if fill == 1
        disp('with filled interior');
        system('./triangulate 1');
    elseif fill == 0
        disp('with unfilled interior');
        system('./triangulate 0');
    elseif fill == 2
        disp('in 2D');
        system('./triangulate 0 2');
    end

    
    %%  Load output
    %   after triangulation, vertices coordinate starts at 0.
    fprintf('Reading Vertices...\n');
    fp = fopen('vert.bin','r');
    vert = fread(fp, [4 inf], 'double')';
    fclose(fp);
    if fill == 2
        tmp = vert(:, [1 2 4]);
        vert = tmp;
        clear tmp;
    end
    
    fprintf('Reading Edges...\n');
    fp = fopen('edge.bin','r');
    edge = fread(fp, [2 inf], 'int32')';
    fclose(fp);
    
    
    fprintf('Processing Edges...\n');
    len = length(vert);    
    if size(edge,1) ~= 0
        tmpgraph = sparse(edge(:,1),edge(:,2),ones(length(edge),1), len, len);
    else
        warning('no edges, skipped');
        skip = 1;
        return;
    end
    tmpgraph = tmpgraph + tmpgraph';
    [I J K] = find(triu(tmpgraph));
    edge = [I J];
    clear I; clear J; clear K;clear tmpgraph;
    
    fprintf('Reading Triangles...\n');
    fp = fopen('triangle.bin','r');
    triangles = fread(fp, [3 inf], 'int32')';
    fclose(fp);
    
    
    % fprintf('Reading Tetrahedrons...\n');
    % fp = fopen('tetrahedron.bin','r');
    % tet = fread(fp, [4 inf], 'int32')';
    % fclose(fp);
    
    
    %%  Write input file for discrete morse
    %   edges and triangles are originally using vertex index starting from 1.
    fprintf('Writing simplices...\n');
    %   since it starts with 0, we put them back to align with original
    %   input.
    vert(:,1) = vert(:,1) + trans(1) + 1;
    vert(:,2) = vert(:,2) + trans(2) + 1;
    vert(:,3) = vert(:,3) + trans(3) + 1;
    write_output(filename,vert, edge-1, triangles-1);
    % write_tetra([filename 'tet_'], vert, tet - 1);
    fprintf('All done!\n');
end

