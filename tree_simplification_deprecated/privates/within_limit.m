function rtn = within_limit(idx,vert)
    lmt = 8*sqrt(2);
    s = vert(idx(1), 1:3);
    t = vert(idx(end), 1:3);
    rtn = true;
    for i = 2:length(idx)-1
        m = vert(idx(i), 1:3);
        std = t - s; norm_std = std ./ vec_len(std);
        var = m - s;
        
        proj_var = (std.* sum(var.*norm_std));
        if vec_len(proj_var - var) > lmt
            rtn = false;
            return;
        end
    end