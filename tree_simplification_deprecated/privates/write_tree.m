function write_tree(tree, vert, path, fname)
    [I J K] = find(triu(tree));
    fpv = fopen([path fname{1}], 'w');
    fpe = fopen([path fname{2}], 'w');
    
    fprintf(fpe, '%d %d %d\n', [I J K]');
    fprintf(fpv, '%d %d %d %f %d\n', vert');
    
    fclose(fpe);
    fclose(fpv);
end