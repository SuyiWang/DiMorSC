'''
DiMorSC Pipeline
'''


from backend import split, DiMorSC
from dataio.struct import Data_Pointer
from dataio import jsonHelper as json
from dataio import fileWriter
import sys

def log_files(filelist):
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
        json_struct = json.parse_all(config)

        if json_struct is not None:
            if json_struct["type"] == "split":
                files = split.pipeline(ptr, json_struct["data"])
                log_files(files)
            elif json_struct["type"] == "DiMorSC":
                files = DiMorSC.pipeline(ptr, json_struct["data"])
                log_files(files)
            else:
                print("[Run]\tInvalid config filetype")
        else:
            print("[Run]\tFailed loading config files")
            