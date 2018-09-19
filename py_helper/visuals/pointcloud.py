import numpy as np
import sys

from vispy import app, visuals, scene, color

def draw(vert, prt=None):
	# tc: visualcontainer
	# if has 4-th dimension, draw with color
	# otherwise, plot grey
	density = vert[:, 3]
	cnorm = density / abs(np.amax(density))
	colors = color.get_colormap("hsl").map(cnorm)#.reshape(density.shape + (-1,))
	V = vert[:, 0:3]

	Scatter3D = scene.visuals.create_visual_node(visuals.MarkersVisual)
	scatter = Scatter3D(parent=prt)
	scatter.set_data(V, edge_color=None, face_color = colors, size=3)
	return scatter


def BuildScene():
	# build canvas
	canvas = scene.SceneCanvas(keys='interactive', title='Graph3D', show=True)


	# Add a ViewBox to let the user zoom/rotate
	view = canvas.central_widget.add_view()
	view.camera = 'turntable'
	view.camera.fov = 45
	view.camera.distance = 1000

	# build visuals
	#GraphPlot = scene.visuals.Line()
	#GraphPlot.parent = view.scene
	return view


def load(filename):
	# load binary simplex file as a graph
	# the file is used for DiMorSC input.
	file = open(filename, 'rb')

	vert_num = np.fromfile(file, dtype=np.int32, count = 1, sep='')
	if vert_num == 0:
		file.close()
		return None, False
	vert =  np.fromfile(file, dtype=np.double, count = vert_num[0] * 4, sep='')
	vert = np.reshape(vert, (-1, 4))
	vert = vert[:, 0:4]

	maxx = int(1e6)
	if len(vert) > maxx:
		print('[point cloud]\tdownsampling input')
		downsample = True

		#get random samples without replacement
		idx = np.random.choice(len(vert), maxx, replace = False)
		# truncate data
		vert = vert[idx, :]
	else:
		downsample = False
	file.close()
	return vert, downsample


if __name__ == '__main__':
	if sys.flags.interactive != 1:
		if len(sys.argv) == 2:
			vert = Load(sys.argv[1])
		else:
			vert = Load('0_dens.bin')
		view = BuildScene()
		draw(vert, view.scene)
		app.run()
