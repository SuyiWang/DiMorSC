'''
[CLASS] Data_Pointer
a data portal used in front end/back end
and _path stores real data pointer for backend.
In front end, Data_Pointer handles only ONE visual object

Properties:
    (numpy.ndarray) _visdata:
        contains data for visualization. Memory consuming
    (string) _type:
        "graph/sc/volume/cloud/json"
        sc is visualized as graph
    (string) _path:
        absolute path of data source
    (int) _ID:
        id for passing information
    (bool) _dataloaded:
        True if the data is completely loaded
    (string) _modtime:
        Last modification time and date of the loaded data
Methods:
    (string) get_full_name():
        get name.ext from path
    (string) get_name():
        get name without extension
    (void) setID():
    (bool) update():
        refetch the data for updates
'''

from dataio import imageLoader
from dataio import jsonHelper as json
import os
import glob

class Data_Pointer():
    # [Types] string:path, numpy.ndarray:visdata, string:objtype, 
    #         int:ID, string:modtime
    # objtype: graph/cloud/volume/sc/json (handled separately)
    def __init__(self, path, visdata={}, objtype=None, offset=[0,0,0]):
        # path should be FULL PATH/abspath
        self._path = path
        self._visdata = visdata
        self._type = objtype
        self._offset = offset

    def _get_name(self):
        return get_name(self._path)
    def _get_dir(self):
        return get_dir(self._path)
    def _get_full_name(self):
        return get_full_name(self._path)
    def _get_ext(self):
        return get_ext(self._path)

    @property
    def type(self):
        return self._type
    @type.setter
    def type(self, value):
        self._type = value

    @property
    def path(self):
        return self._path
    @path.setter
    def path(self, value):
        if value[-1] == '/':
            self._path = value[0:-1]
        else:
            self._path = value

    @property
    def offset(self):
        return self._offset
    @offset.setter
    def offset(self, value):
        self._offset = value


def get_name(path):
    fullname = get_full_name(path)
    namesplit = fullname.split('.')
    if len(namesplit) > 1:
        return namesplit[-2]
    else:
        return namesplit[-1]

def get_dir(path):
    return os.path.dirname(path)

def get_full_name(path):
    return path.split('/')[-1]

def get_ext(path):
    fullname = get_full_name(path)
    namesplit = fullname.split('.')
    if len(namesplit) <= 1:
        return None
    else:
        return namesplit[-1]
