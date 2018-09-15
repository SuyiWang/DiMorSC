'''
DiMorSC Pipeline

Methods:
    public methods:
    (list) pipeline(Data_Pointer data_ptr, json-like parameters, string workpath)
    (string) merge(Data_Pointer data_ptr, string workpath)
    
    private methods:
    (string) _Gsmooth(numpy.ndimage data, double sigma, double thd,
    string filename):
        smoothes data and write result to filename

    (string) _triangulate(int id, int fill, int dim):
        triangulate point cloud and write result as simplicial complex

    (string) _run(int id, double persist, int dim):
        Run DiMorSC on output/id.sc and simplify according to 'persist'
'''


import subprocess
from dataio import fileWriter, imageLoader, fileReader
from dataio.struct import Data_Pointer
from dataio import jsonHelper as json
import scipy.ndimage as ndimage
import os


# set to None to show subprocess output
# set to subprocess.DEVNULL to hide them
stdoutREDIRECT = None#subprocess.DEVNULL


def pipeline(data_ptr, parameters, workpath="output/"):
    '''
    Process a single file through pipeline:
        smooth, triangulate, trace, tree_simplify
    '''

    # if workfolder does not exist, create it
    if workpath[-1] != '/':
        workpath = workpath + '/'
    os.makedirs(workpath, exist_ok=True)

    # maintains a data_ptr to be processed
    now_ptr = data_ptr
    # list of generated files
    rtn = []

    for para in parameters:
        if para['action'] == 'preprocess':
            #prepare
            print('    [DiMorSC]\tLoading data:')
            data = imageLoader.load(now_ptr.path)
            
            sigma = float(para['sigma'])
            cut_thd = float(para['threshold'])
            print('    [DiMorSC]\tSmoothing:sigma %.2f, threshold %.2f'
                  %(sigma,cut_thd))
            
            _check_offset(data.shape, now_ptr)

            # run
            file = _Gsmooth(
                        data, 
                        sigma, 
                        cut_thd,
                        workpath + now_ptr._get_name(),
                        now_ptr.offset[0:3]
                        )
            
            # collect
            rtn.append(('dens', file))
            now_ptr = Data_Pointer(file, objtype='cloud')
            del data

        if para['action'] == 'triangulation':
            # prepare
            print('    [DiMorSC]\ttriangulating target:', now_ptr.path)

            # run
            file = _triangulate(now_ptr.path)

            # collect
            now_ptr = Data_Pointer(file, objtype='sc')
            rtn.append(('sc', file))

        if para['action'] == 'DiMorSC':
            # prepare
            print('    [DiMorSC]\trunning DiMorSC on:', now_ptr.path)
            pthd = float(para['threshold'])
            print('        Persistent Threshold:%f'%pthd)

            # run
            file = _run(now_ptr.path, pthd)

            # collect
            name, ext = os.path.splitext(os.path.basename(file))
            graphcfg = [
            ''.join([name, '_vert.txt']),
            ''.join([name, '_edge.txt']),
            ]
            ininame = ''.join([file, '.ini'])
            fileWriter.write_string_list(graphcfg, ininame)
            rtn.append(('ini', ininame))
            now_ptr = Data_Pointer(ininame, objtype='graph')

        if para['action'] == 'to_tree':
            # prepare input config, append configuration lines
            print('    [DiMorSC]\tgraph2tree')
            _append_config(now_ptr.path, para)
            
            # run
            files = _graph2tree(now_ptr.path)

            # collect
            for file in files:
                rtn.append(('graph', file))

    print('    files to be collected:')
    filelist = '\n'.join([('\t'+x[0]+': '+x[1]) for x in rtn])
    print(filelist)
    return rtn


def merge(data_ptr, workpath='output/'):
    # assuming when merge, 'pos' in 'data' must have full block info.
    # this should be checked when processing each block in previous stage.
    if workpath[-1] != '/':
        workpath = workpath + '/'
    os.makedirs(workpath, exist_ok=True)
    file_structure = json.parse_all(data_ptr.path, jsontype = 'block')
    
    # write input for c++ merger
    filename = _write_merger_config(file_structure['data'], 
                                   file_structure['overlap'], 
                                   workpath)
    # call merger
    subprocess.call(["./bin/merge_graph", filename],
                    stdout=stdoutREDIRECT)
    return [('sc', "merged"),]


def _rmvExt(filename):
    return filename.rsplit('.', 1)[0]


def _Gsmooth(data, sigma = 2, thd = 0.01, filename = 'output/0', trans = (0,0,0)):
    smooth_data = ndimage.filters.gaussian_filter(data, sigma)
    #data = ndimage.filters.minimum_filter(data, [5,5,3])
    #data = ndimage.filters.uniform_filter(data, [7,7,3])
    
    fileWriter.write_bin(smooth_data, filename, thd, trans)
    return filename + '.dens'


def _triangulate(filename='output/0.dens', fill=0, dim=3):
    subprocess.call(["./bin/Triangulate",
          filename,
          str(fill),
          str(dim)
         ],
         stdout=stdoutREDIRECT)
    return _rmvExt(filename) + '.sc'


def _run(filename='output/0.sc', persist=1, dim=3):
    outputprefix = _rmvExt(filename)
    subprocess.call(["./bin/DiMorSC", 
          filename,
          outputprefix,
          str(persist),
          str(dim)
         ],
         stdout=stdoutREDIRECT)
    return outputprefix
    # 1e7 points cost about 32G memory and 15min time.


def _graph2tree(ininame = 'output/0.ini'):
    rtn = []
    subprocess.call(["./bin/graph2tree", ininame])
    # for now it did not return anything
    # the result can be viewed using view.py
    return rtn


def _write_merger_config(data, overlap, workpath):
    filename = workpath + 'merger_config'
    with open(filename, "w") as f:
        # output filename
        f.write("merged\n")
        # overlap info
        overlap_str = [str(x) for x in overlap]
        f.write(" ".join(overlap_str))
        f.write("\n");
        for item in data:
            # block file prefix, this is usually integer
            f.write(item["file"])
            f.write(" ")
            
            # 1. block locator - converted to XYZ order
            # This offset should be processed in reversed order 
            # to transit from ZYX system
            # 2. convert the positions to int - In current version
            # offset is based on pixels, no need of float
            pos_str = [str(int(x)) for x in reversed(item["pos"][0:3])]
            pos_str += [str(int(x)) for x in reversed(item["pos"][3:6])]
            f.write(" ".join(pos_str))
            f.write("\n")
    return filename


def _check_offset(shape, obj):
    # inner coordinate, offset are zyx order
    if len(obj.offset) == 6:
        pass
    elif len(obj.offset) == 3:
        ending = [obj.offset[0]+shape[0],
                  obj.offset[1]+shape[1],
                  obj.offset[2]+shape[2]]
        obj.offset += ending
    else:
        print("    [DiMorSC]\tIncorrect offset. offset must be 3/6 integers")
    return


def _append_config(inifile, para_tree):
    # get first two lines from ini file
    lines = fileReader.parse_ini(inifile, line_count=2)
    name, ext = os.path.splitext(inifile)

    lines.append(name)

    root = para_tree.get("root", "0 0 0")
    lines.append(root)

    saddle = para_tree.get("saddle", 0)
    lines.append(str(saddle))

    component = para_tree.get("component", 0)
    lines.append(str(component))

    fileWriter.write_string_list(lines, inifile)



def write_files(filelist):
    fileWriter.write_string_list([x[1] for x in filelist],
                                 'output/log.txt')


if __name__  == "__main__":
    if len(sys.argv) != 3:
        print ("[Run]\tusage: run <path_to_data> <json_file>")
        print ("See data/parameters.json for an example")
    else:
        target = sys.argv[1]
        config = sys.argv[2]

        # validate target
        ptr = Data_Pointer(target)
        # validate config
        para = json.parse_data(config, jsontype = 'DiMorSC')

        if para is not None:
            files = process(ptr, para)
            # write all files to 'output/log.txt'
            write_files(files)
