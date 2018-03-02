function [skip] = Triangulate( filename, fill, trans, id )
    skip = 0;
    filename = [filename int2str(id)];
    
    %%  Triangulate - 0) do not fill the cube; 1) fill the cube interier with a tetrahedron
%     system('rm vert.bin');
%     system('rm edge.bin');
%     system('rm triangle.bin');
    fprintf('Trianguating...');
    if fill == 1
        disp('with filled interior');
        system(['./core/triangulate ' int2str(id) ' 1 3']);
    elseif fill == 0
        disp('with unfilled interior');
        system(['./core/triangulate ' int2str(id) ' 0 3']);
    elseif fill == 2
        disp('in 2D');
        system(['./core/triangulate ' int2str(id) ' 0 2']);
    end
	
	delete([int2str(id) '_dens.bin']);
    
    %%  Load output
    %   after triangulation, vertices coordinate starts at 0.
    fprintf('Reading Vertices...\n');
    fp = fopen([int2str(id) '_vert.bin'],'r');
    vert = fread(fp, [4 inf], 'double')';
    fclose(fp);
	delete([int2str(id) '_vert.bin']);
    if fill == 2 && length(vert) > 0
        tmp = vert(:, [1 2 4]);
        vert = tmp;
        clear tmp;
    end
    
    fprintf('Reading Edges...\n');
    fp = fopen([int2str(id) '_edge.bin'],'r');
    edge = fread(fp, [2 inf], 'int32')';
    fclose(fp);
    delete([int2str(id) '_edge.bin']);
    
    fprintf('Processing Edges...\n');
    len = length(vert);    
    if size(edge,1) ~= 0
        tmpgraph = sparse(edge(:,1),edge(:,2),ones(length(edge),1), len, len);
    else
        warning('no edges, skipped');
        skip = 1;
        write_output(filename);
		delete([int2str(id) '_triangle.bin']);
        return;
    end
    tmpgraph = tmpgraph + tmpgraph';
    [I J K] = find(triu(tmpgraph));
    edge = [I J];
    clear I; clear J; clear K;clear tmpgraph;
    
    fprintf('Reading Triangles...\n');
    fp = fopen([int2str(id) '_triangle.bin'],'r');
    triangles = fread(fp, [3 inf], 'int32')';
    fclose(fp);
	delete([int2str(id) '_triangle.bin']);
    
    
    % fprintf('Reading Tetrahedrons...\n');
    % fp = fopen('tetrahedron.bin','r');
    % tet = fread(fp, [4 inf], 'int32')';
    % fclose(fp);
    
    
    %%  Write input file for discrete morse
    %   edges and triangles are originally using vertex index starting from 1.
    fprintf('Writing simplices...\n');
    %   since it starts with 0, we put them back to align with original
    %   input.
    vert(:,1) = vert(:,1) + double(trans(1) + 1);
    vert(:,2) = vert(:,2) + double(trans(2) + 1);
    vert(:,3) = vert(:,3) + double(trans(3) + 1);
    write_output(filename,vert, edge-1, triangles-1);
    % write_tetra([filename 'tet_'], vert, tet - 1);
    fprintf('All done!\n');
end

