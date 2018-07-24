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

# Read volume
def Load(filename):
    print ("[visuals.volume] \tvolume.load TODO")

def Draw(data, tc):
	# tc: visualcontainer
	#   More reference: vispy/examples/basics/visuals/image_transforms.py
	
	# translate group 1
	# vol1 = np.rot90(data, 1, axes=(2,1))
	
	# translate group 2
	vol1 = np.rot90(data, 1, axes = (2, 0))
	vol1 = np.flip(vol1, axis = 2)

	shp = vol1.shape
	if shp[2] % 2 == 1:
		warnings.warn('[volume] \tvolume dimension 2 has odd number of size, trying to fix', RuntimeWarning)
		zo = np.zeros((shp[0], shp[1], 1), dtype='f8')
		vol1 = np.concatenate((vol1, zo), axis = 2)
	#vol1 = data[:,::-1,:]
	#vol1 = vol1[:,:,::-1]
	
	# Set whether we are emulating a 3D texture
	emulate_texture = False
	# Create the volume visuals, only one is visible
	volume1 = scene.visuals.Volume(vol1, parent=tc.view.scene, threshold=0.001,
	                               emulate_texture=emulate_texture)
	#volume1.transform = scene.STTransform(translate=(64, 64, 0))
	#transform = scene.MatrixTransform()
	#transform.rotate(270, (0,0,1))
	#volume1.transform = transform
	return volume1

# volume has method set_data, which can be used to update data.
def new(obj, viewer):
	volume = Draw(obj._visdata['volume'], viewer)
	viewer.visuallist[obj._ID] = volume
def update(obj, viewer):
	viewer.visuallist[obj._ID].set_data(obj._visdata['volume'])
