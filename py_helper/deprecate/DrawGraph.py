import numpy as np
import sys

from vispy import app, visuals, scene

def LoadGraph(fileprefix):

	vertinput = np.loadtxt(fileprefix+'_vert.txt')
	vert = vertinput[:, 0:3]
	edgeinput = np.loadtxt(fileprefix + '_edge.txt')
	edge = edgeinput[:, 0:2].astype('i')
	# now vertex are indexed from 1. Hence we convert it back.
	edge = edge-1
	
	return edge, vert

def DrawGraph(edge, vert, targetscene):
	scene.visuals.Line(pos=vert, color='#f006', method='gl', connect=edge, parent=targetscene.scene)
	return


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


if __name__ == '__main__':
    if sys.flags.interactive != 1:
    	edge, vert = LoadGraph('0')
    	view = BuildScene()
    	Draw(edge, vert, view.scene)
        app.run()