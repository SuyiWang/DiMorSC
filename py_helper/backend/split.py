'''
This pipeline is for 
    splitting neurons in very large image into smaller blocks
    tracing neurons contained in multiple blocks
    merge splitted neurons

Methods:
    (list[string]) pipeline(Data_Pointer obj, para parameters):
        process obj with 'parameters' extracted from parameter box
'''

from backend import DiMorSC
from dataio.struct import Data_Pointer
from dataio import jsonHelper as json
from dataio import splitter

# detect involved files into Data_Pointer
def detect_files(data_ptr):
    if data_ptr.type!='json':
        # if contains single file, pass data_ptr as a list
        return [data_ptr, ]
    else:
        # if not, create fake Data_Pointer for each file
        file_structure = json.parse_data(data_ptr.path, jsontype = 'block')
        rtn = []
        path = data_ptr._get_dir()
        for block in file_structure:
            rtn.append(Data_Pointer(path+'/'+block['file'],
                               offset=block['pos'])
                      )
        return rtn

# process Data_Pointer passed from GUI
def pipeline(data_ptr, parameters, workpath='output/'):
    if workpath[-1] != '/':
        workpath = workpath + '/'
    # list of objects to be visualized
    rtn = []
    # maintains a Data_Pointer to be processed
    now_ptr = data_ptr
    data_ptr_list = None

    # batch process to be handled
    for para in parameters:
        if para['action'] == 'split':
            print('[Backend]\t spliting target:', now_ptr.path)
            print('    split size:\t', para['size'])
            # a folder "[objname]_split" will be created in work path
            # practically, output file is "[objname].json" in that folder
            file = splitter.to_blocks(now_ptr, 
                                      para['size'], 
                                      workpath=workpath)

            # create a new Data_Pointer
            now_ptr = Data_Pointer(file,
                                  objtype='json')
            rtn.append(('json', file))

        if para['action'] == 'trace':
            print('[Backend]\tTrace target:', now_ptr.path)
            # Read file list to be processed
            # returns a list of files, len(list) may be 1
            data_ptr_list = detect_files(now_ptr)
            # a folder with dataname will be created in workpath
            # to contain all outputs
            OFFSET_UPDATE = False
            dataname = now_ptr._get_name()
            for item in data_ptr_list:
                if not OFFSET_UPDATE and len(item.offset) != 6:
                    OFFSET_UPDATE = True
                tracefile = DiMorSC.pipeline(item, 
                                        para['children'], 
                                        workpath+dataname)

            # rewrite if json file does not contain full block info
            # only update for json
            if OFFSET_UPDATE:
                _update_offset(now_ptr, data_ptr_list)
            # may allow collect intermediate result

            # if file is collected here, ignore it in the future
            if len(data_ptr_list) == 1:
                rtn = rtn + tracefile
                now_ptr = Data_Pointer(tracefile, objtype='oldgraph')
                tracefile = None

        # collect results
        if para['action'] == 'merge':
            # if object is json,
            # the function will search in "workpath+dataname"
            # for the tracing result of blocks.
            # if object is a file/folde, it creates a folder in workpath
            # AVOID putting data file in workpath
            if now_ptr.type == 'json':
                # merge fragments if contains multiple blocks
                print('[Backend]\tMerging target:', now_ptr.path)
                dataname = now_ptr._get_name()
                tracefile = DiMorSC.merge(now_ptr, workpath+dataname)
                now_ptr = Data_Pointer(tracefile, objtype='oldgraph')
                rtn.append(('oldgraph', tracefile))
            else:
                # do not change now_ptrect
                print('[Backend]\tCollecting single graph for target:',
                      now_ptr.path)
                if 'tracefile' in locals() and tracefile is not None: 
                    now_ptr = Data_Pointer(tracefile, objtype='oldgraph')
                    rtn.append(('oldgraph', tracefile))

        if para['action'] == 'to_tree':
            print('[DiMorSC]\tgraph2tree')
            # prepare input config

            files = _graph2tree(now_ptr.path)
            for file in files:
                rtn.append(('graph', file))

    # return all files
    return rtn


def _update_offset(obj, data_ptr_list):
    if obj.type == 'json':
        print('[Backend]\tUpdating offset in:', obj._path)
        jsonoutput = json.parse_all(obj.path, jsontype = 'block')
        file_structure = jsonoutput['data']
        assert len(file_structure) == len(data_ptr_list), "[Backend]\tlength of blocks mismatch"
        for i in range(len(data_ptr_list)):
            file_structure[i]['pos'] = data_ptr_list[i].offset
        
        json.write(jsonoutput, obj._path)


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
        para = json.parse_data(config, jsontype = 'split')

        if para is not None:
            files = pipeline(ptr, para)
            # write all files to 'output/log.txt'
            write_files(files)
