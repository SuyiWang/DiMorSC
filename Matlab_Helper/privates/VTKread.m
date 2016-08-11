function [field,value] = VTKread(filename, verbose)

% read_vtk - read data from VTK file.
%
%   [vertex,face] = read_vtk(filename, verbose);
%
%   'vertex' is a 'nb.vert x 3' array specifying the position of the vertices.
%   'face' is a 'nb.face x 3' array specifying the connectivity of the mesh.
%
%   Copyright (c) Mario Richtsfeld

if nargin<2
    verbose = 1;
end

fid = fopen(filename,'r');
if( fid==-1 )
    error('Can''t open the file.');
    return;
end

str = fgets(fid);   % -1 if eof
if ~strcmp(str(3:5), 'vtk')
    error('The file is not a valid VTK one.');    
end

%%% read header %%%
str = fgets(fid);
str = fgets(fid);
str = fgets(fid);
str = fgets(fid);
fieldcount = 3;
while str ~= -1
    NAME = sscanf(str, '%s', 1);
    if strcmp(NAME,'POINTS')
        % read vertices
        nvert = sscanf(str,'%*s %d', 1);
        [A,cnt] = fscanf(fid,'%f %f %f', 3*nvert);
        if cnt~=3*nvert
            warning('Problem in reading vertices.');
        end
        A = reshape(A, 3, cnt/3);
        points = A;
        str = fgets(fid);
    elseif strcmp(NAME,'LINES')
        nvert = sscanf(str,'%*s %d %d', 2);
        I = zeros(nvert(2)-nvert(1)*2,1);
        J = zeros(nvert(2)-nvert(1)*2,1);
        tot = 1;
        
        % assume cnt>=3
        for i = 1:nvert(1)
            str = fgets(fid);
            [A,cnt] = sscanf(str,'%d');
            I(tot:tot+cnt-3) = A(2:cnt-1);
            J(tot:tot+cnt-3) = A(3:cnt);
            tot = tot + cnt - 2;
        end
        
    elseif strcmp(NAME,'VERTICES')
        nvert = sscanf(str,'%*s %d', 1);
        [A,cnt] = fscanf(fid,'%*f %f', nvert);
        if cnt~=nvert
            warning('Problem in reading vertices.');
        end
        vertices = A;
        str = fgets(fid);
    elseif strcmp(NAME,'POINT_DATA')
        str = fgets(fid);
        str = fgets(fid);
    	continue;
    elseif strcmp(NAME,'CELL_DATA')
        break;
    else
        field{fieldcount} = NAME;
        nvert = sscanf(str,'%*s %*d %d', 1);
        [A,cnt] = fscanf(fid,'%f', nvert);
        if cnt~=nvert
            warning('Problem in reading data field.');
        end
        value{fieldcount} = A;
        fieldcount = fieldcount+1;
        str = fgets(fid);
    end
    
    
    str = fgets(fid);
end
fclose(fid);
G = sparse(I+1, J+1, ones(length(I),1), length(points), length(points));

[G, vert, mergemap, idx] = mergeduplicate(G, points');

for i = 3:9
    value{i} = value{i}(idx);
end
% vertices = vertices+1;
field{1} = 'G'; field{2} = 'vert';
value{1} = G; value{2} = vert;
