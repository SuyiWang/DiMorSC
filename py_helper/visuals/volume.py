from itertools import cycle
import warnings
import numpy as np

from vispy import app, io
from vispy import scene

from vispy.color import get_colormaps, BaseColormap
from vispy.visuals.transforms import STTransform, MatrixTransform

# Reference: http://vispy.org/visuals.html
# volume.set_data will result in a COPY of the volume data.
# In general, we need to keep both, 
# because the copied data will be changed to a form for rendering

def draw(data, prt=None, translate=[0,0,0]):
    # sometimes the volume coordinates are inconsistent with the actual one
    # see binWriter.py for currently found correct one.
    # if not, manually adjust with the following traslate groups.

    shp = data.shape
    if shp[2] % 2 == 1:
        warnings.warn('[volume] \tvolume dimension 2 has odd number of size, trying to fix', RuntimeWarning)
        zo = np.zeros((shp[0], shp[1], 1), dtype='f8')
        data = np.concatenate((data, zo), axis = 2)
    
    # Set whether we are emulating a 3D texture
    emulate_texture = False
    # Create the volume visuals, only one is visible
    volume1 = scene.visuals.Volume(data, parent=prt,
                                   threshold=0.001,
                                   emulate_texture=emulate_texture,
                                   method='mip')
    
    # translate coordinate is provided in ZYX
    # but screen coordinate is in XYZ, so we reverse it
    translate = translate[::-1]
    volume1.transform = scene.STTransform(translate=translate)
    return volume1

def multi_draw(blocklist, prt=None):
    # obj path is json. visdata is a list of vobject of type volume.
    # get data.
    
    # draw for each item in blocklist a visual
    volume_block_visual = []
    for vobj in blocklist:
        vol_visual = draw(vobj._visdata['volume'],
                          prt=prt,
                          translate=vobj.offset)
        volume_block_visual.append(vol_visual)

    # data list is a list of 2-tuple (data, offset)
    return  volume_block_visual
