'''
DiMorSC Pipeline

Methods:
	(string) Gsmooth(numpy.ndimage data, double sigma, double thd,
	string filename):
		smoothes data and write result to filename

	(string) triangulate(int id, int fill, int dim):
		triangulate point cloud and write result as simplicial complex

	(string) run(int id, double persist, int dim):
		Run DiMorSC on output/id.sc and simplify according to 'persist'
'''


from subprocess import call
from dataio import binWriter
import scipy.ndimage as ndimage
from dataio import imageLoader
import os
from dataio import config
from glob import glob


def rmvExt(filename):
	return filename.rsplit('.', 1)[0]


def loadreal(data_pointer):
	if os.path.isdir(data_pointer):
		data = imageLoader.Load(data_pointer, 'tif')
	else:
		data = imageLoader.Load(data_pointer)
	return data

def Gsmooth(data, sigma = 2, thd = 0.01, filename = 'output/0'):
	smooth_data = ndimage.filters.gaussian_filter(data, sigma)
	#data = ndimage.filters.minimum_filter(data, [5,5,3])
	#data = ndimage.filters.uniform_filter(data, [7,7,3])
	binWriter.write(smooth_data, filename, thd)
	return filename + '.dens'


def triangulate(filename='output/0.dens', fill=0, dim=3):
	# writes to output/
	call(["./bin/Triangulate", filename, str(fill), str(dim)])
	return rmvExt(filename) + '.sc'


def DiMorSC(filename='output/0.sc', persist=1, dim=3):
	# writes to path given by output prefix
	# here it's the same directory as input 
	outputprefix = rmvExt(filename)
	call(["./bin/DiMorSC", filename, outputprefix, str(persist), str(dim)])
	return outputprefix
	# 1e7 points cost about 32G memory and 15min time.

def graph2tree(ininame = 'output/0.ini'):
	rtn = []
	inicfg = config.getGraphConfig(ininame)
	if inicfg is not None:
		outputprefix = ''.join([inicfg[2], '_tree'])
		call(["./bin/graph2tree", ininame])
	# for now it did not return anything
	# the result can be viewed using view.py
	return rtn
