

def parse_ini(filename, line_count=None):
    rtn = []
    with open(filename) as fp: 
        for line in fp:
            if line_count is not None and line_count <= 0:
                break
            
            # should consider comments later
            line_count -= 1

            # remove \n from string
            rtn.append(line[0:-1])
    return rtn
