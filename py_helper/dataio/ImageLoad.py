#OLD VERSION
# python -W ignore ImageLoad.py Olfactory_OP7_trunc.tif --tif

import numpy as np
import scipy.ndimage as ndimage
from sys import argv
import glob, os
from skimage import io
from PIL import Image


import re
def tryint(s):
    try:
        return int(s)
    except:
        return s
 
def alphanum_key(s):
    """ Turn a string into a list of string and number chunks.
        "z23a" -> ["z", 23, "a"]
    """
    return [ tryint(c) for c in re.split('([0-9]+)', s) ]


class ImageLoader():
	
	def Load(self, filepath, suffix=''):
		if os.path.isdir(filepath):
			if filepath[-1] is not '/':
				filepath = filepath + '/'
			files = self.GetFilelist(filepath, suffix)
			data = []
			for file in files:
				data.append(self.LoadSingleTIF(filepath+file))

			rtn = np.stack(data, axis=2)
			print ""
			print "Files stacked with size"
			print rtn.shape
			return rtn

		elif os.path.isfile(filepath):
			#  This should be changed to consider more file types
			return self.LoadTIF(filepath)
		else:
			print "Unable to load %s" %filepath
			return None

	def GetFilelist(self, path, suffix):
		print "Reading files in "+path
		files = os.listdir(path)
		return sorted(files, key=alphanum_key)


	def LoadTIF(self, filename):
		#files = GetFilelist(path, "tif")
		'''
		# for loading a single simage
		im = Image.open(filename)
		imarray = numpy.array(im)
		print imarray.shape
		print im.size
		'''
		
		# for loading image stacks
		printfilename = filename.split('/')[-1]
		if filename is None:
			print "LoadTIF: Nothing to do."
			return None
		im  = io.imread(filename)
		if im is None:
			print "Failed to load %s"%filename
			return None
		#  [TODO?] coule be changed for better precision
		#  F8 is 64bit = double

		if im.dtype != np.uint8:
			maxx = np.amax(im)
			#print "max value: " + str(maxx)
			rtn = np.array(im, dtype='f8')
			rtn = rtn*(255.0/maxx)
			###Vispy has no constraints for uint8 or uint16

		else:
			rtn = np.array(im, dtype='f8')
		
		msg = "Image %s loaded. Size: "%printfilename  + str(rtn.shape) + '\r'
		print msg ,
		return rtn

	def LoadSingleTIF(self, filename):
		#files = GetFilelist(path, "tif")
		'''
		# for loading a single simage
		im = Image.open(filename)
		imarray = numpy.array(im)
		print imarray.shape
		print im.size
		'''
		
		# for loading image stacks
		if filename is None:
			print "LoadTIF: Nothing to do."
			return None
		im  = io.imread(filename)
		if im is None:
			print "Failed to load %s"%filename
			return None
		#  [TODO?] coule be changed for better precision
		#  F8 is 64bit = double
		
		if (im.ndim > 2):
			im = np.sum(im, axis=2)
			print "ndim > 2"
		if im.dtype != np.uint8:
			maxx = np.amax(im)
			input = np.array(im, dtype='f8')
			input = input*(255.0/maxx)
			print np.amax(input), maxx
			print "not uint8"
		else:
			input = np.array(im, dtype='f8')
		msg = "Image %s loaded. Size: "%filename  + str(input.shape) + '\r'
		print msg ,
		return input

	def imagetestwrite(self, A, filename):
		A = np.array(A, dtype=np.uint8)
		im = Image.fromarray(A)
		file = filename.split('/')[-1]
		file = file.split('.')[0]
		im.save('test/'+file+'.png','PNG')

	def getprefix(filename):
		return 

if __name__  == "__main__":
	print "parameters detected: " + str(len(argv))
	if argv[2] == '--tif':
		IL = ImageLoader()
		data = IL.Load(argv[1], 'tif')
	else:
		print "usage: ImageLoad <filename> <filetype>"

		