import numpy as np
import sys
import os

from vispy import app, visuals, scene, color

# Load file with no extension - as graph
def Load(fileprefix):
	vertinput = np.loadtxt(fileprefix+'_vert.txt')
	if len(vertinput) > 0:
		vert = vertinput[:, 0:4]
	else:
		vert = None

	edgeinput = np.loadtxt(fileprefix + '_edge.txt')
	if len(edgeinput) > 0:
		edge = edgeinput[:, 0:2].astype('i')
		# now vertex are indexed from 1. Hence we convert it back.
		edge = edge-1
	else:
		edge = None
	
	return edge, vert

# Done
def loadSWC(filename):
	swcinput = np.loadtxt(filename)
	if len(swcinput) > 0:
		vert = swcinput[:, 2:6]
		# assuming the first line is root.
		edge = swcinput[1:, (0, 6)].astype('i')
		edge = edge - 1
	else:
		edge = None
		vert = None
	return edge, vert

# Done
def loadini(filename):
	graphname = []
	counter = 0
	with open(filename, "r") as fp:
		line = fp.readline()
		while line:
			if line.startswith('#'):
				line = fp.readline()
				continue
			# readline has \n in the end
			graphname.append(line.rstrip('\n'))
			counter += 1
			if counter == 2:
				break
			line = fp.readline()

	fp.close()

	if len(graphname) < 2:
		print('[graph] not enough filename spceified. doing nothing')
		return None, None

	path = getPath(filename)
	print ('[graph]: Loading from path ' + path)
	vertinput = np.loadtxt(path + graphname[0])
	if len(vertinput) > 0:
		vert = vertinput[:, 0:4]
	else:
		vert = None
	edgeinput = np.loadtxt(path + graphname[1])
	if len(edgeinput) > 0:
		edge = edgeinput[:, 0:2].astype('i')
		# now vertex are indexed from 1. Hence we convert it back.
		edge = edge-1
	else:
		edge = None
	
	return edge, vert


def getPath(s):
	paths = s.split('/');
	if len(paths) == 1:
		return ""
	else:
		return s.rsplit('/', 1)[0] + '/'


def DrawPoints(vert, prt=None):
	# tc: visualcontainer
	# if has 4-th dimension, draw with color
	# otherwise, plot grey
	face_color = (1, 1, 1, .5)

	if vert.shape[1] < 4:
		density = np.ones(vert.shape[0])
	else:
		density = vert[:, 3]
	cnorm = density / abs(np.amax(density))
	colors = color.get_colormap("hsl").map(cnorm)#.reshape(density.shape + (-1,))
	V = vert[:, 0:3]

	Scatter3D = scene.visuals.create_visual_node(visuals.MarkersVisual)
	scatter = Scatter3D(parent=tc.view.scene)
	scatter.set_data(V, edge_color=None, face_color = colors, size=3)
	return scatter



def draw(edge, vert, prt=None):
	# tc: visualcontainer
	vert = vert[:, 0:3]
	graph = scene.visuals.Line(pos=vert,
							   color='#f006',
							   method='gl',
							   connect=edge,
							   parent=prt,
							   width = 3)
	return graph


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


def LoadSimplexGraph(filename):
	# load binary simplex file as a graph
	# the file is used for DiMorSC input.
	file = open(filename, 'rb')

	vert_num = np.fromfile(file, dtype=np.int32, count = 1, sep='')
	if vert_num == 0:
		file.close()
		return None, None
	vert =  np.fromfile(file, dtype=np.double, count = vert_num[0] * 4, sep='')
	vert = np.reshape(vert, (-1, 4))
	vert = vert[:, 0:3]

	edge_num = np.fromfile(file, dtype=np.int32, count = 1, sep='')
	if edge_num == 0:
		file.close()
		return None, None
	edge = np.fromfile(file, dtype=np.int32, count = edge_num[0] * 2, sep='')
	edge = np.reshape(edge, (-1, 2))
	file.close()
	return edge, vert


if __name__ == '__main__':
	if sys.flags.interactive != 1:
		if len(sys.argv) > 1:
			edge, vert = Load(sys.argv[1])
		else:
			edge, vert = Load('0')
		
		view = BuildScene()
		# create fake scene with view attribute
		tc = lambda: None
		setattr(tc, 'view', view)

		DrawPoints(vert, prt=tc.view.scene)
		draw(edge, vert, prt=tc.view.scene)
		app.run()
