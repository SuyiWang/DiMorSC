
from vispy import app
from visuals.VisualManager import MyCanvas
import sys
from dataio.struct import Data_Pointer
from dataio import config

if __name__  == "__main__":
	print ("[viewer]\tparameters detected: " + str(len(sys.argv)))
	if len(sys.argv) == 2:
		vc = MyCanvas()
		vc.setup()
		# load file
		obj = Data_Pointer(sys.argv[1])
		vc.insert(obj)
		app.run()
	elif len(sys.argv) == 3:
		vc = MyCanvas()
		vc.setup()
		# load file list
		files = config.getLines(sys.argv[1])
		for file in files:
			if len(file) < 1:
				continue
			obj = Data_Pointer(file)
			vc.insert(obj)
		app.run()
	else:
		print ("[Run]\tusage: run <data_file> [list]")
