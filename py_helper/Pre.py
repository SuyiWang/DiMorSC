from dataio.ImageLoad import ImageLoader
from dataio.ImageWrite import ImageWriter
import scipy.ndimage as ndimage
from sys import argv
from Exec import Triangulate, DiMorSC
import visuals


#python Pymage/Pre.py Olfactory_Projection_Fibers/Image_Stacks/OP_1/ --render

if __name__  == "__main__":
	print "parameters detected: " + str(len(argv))
	if argv[2] == '--tif':
		IL = ImageLoader()
		data = IL.Load(argv[1], 'tif')

		# canvas, view = visuals.create_canvas()
		# visuals.volume.Draw(data, view.scene)
		# visuals.Show(canvas)
		
		smooth_data = ndimage.filters.gaussian_filter(data, 2)

		ImageWriter(smooth_data, '0', 10)
		Triangulate(0, 0)
		DiMorSC(0, 500)
		
	elif argv[2] == '--render':
		IL = ImageLoader()
		data = IL.Load(argv[1], 'tif')
		canvas, view = visuals.create_canvas()
		visuals.volume.Draw(data, view.scene)
		edge, vert = visuals.graph.Load('0')
		visuals.graph.Draw(edge, vert, view.scene)
		visuals.Show(canvas)
	else:
		print "usage: Pre <filename> <filetype>"