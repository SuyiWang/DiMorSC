'''
DiMorSC Pipeline
'''
from backend import pipeline
import sys
from os.path import basename
from dataio import config

def getStage(filetype):
	if filetype == 'volume':
		return 0
	elif filetype == 'cloud':
		return 1
	elif filename == 'sc':
		return 2
	elif filename == 'graph':
		return 3
	else:
		return -1

def writeFiles(filelist):
	f = open('output/log.txt', 'w')
	if f is not None:
		for item in filelist:
			f.write(item[1])
		f.close()
	else:
		print('[run]\tError creating log to output/log.txt')

'''
0,{inputname}
1,{type}
2,{sigma}
3,{densThd}
4,{persist}
5,{root}
6,{saddleThd}
'''
def process(cfg):
	rtn = []

	# [TODO] This could be changable through parameter box
	workpath = 'output/'
	# decides which input to load
	inputpointer = cfg[0]
	name = basename(cfg[0])
	stage = getStage(cfg[1])

	if stage < 0:
		print("[Run]\tinvalid filetype. Filetype should be volumd/cloud/sc/graph.")
		return []

	if stage < 1:
		# read image data
		data = pipeline.loadreal(inputpointer)
		print('[Run]\tSmoothing:')
		# smooth data and write to point cloud
		file = pipeline.Gsmooth(
					data, 
					float(cfg[2]), 
					float(cfg[3]),
					workpath + name)
		rtn.append(('dens', file))
		inputpointer = file
		del data

	if stage < 2:
		print('[Run]\ttriangulating:')
		file = pipeline.triangulate(inputpointer)
		inputpointer = file
		# usually do not visualize sc, since its big size and 
		# is generally no more useful than the point cloud
		# rtn.append(('sc', file))

	if stage < 3:
		print('[Run]\trunning DiMorSC:')
		file = pipeline.DiMorSC(inputpointer, float(cfg[4]))
		graphcfg = [
		''.join([file, '_vert.txt']),
		''.join([file, '_edge.txt']),
		file,
		cfg[5],
		cfg[6]
		]
		file = ''.join([file, '.ini'])
		config.writeGraphConfig(file, graphcfg)
		rtn.append(('ini', file))
		inputpointer = file

	if stage < 4:
		print('[Run]\tgraph2tree')
		# if file does not exist, create an example
		files = pipeline.graph2tree(inputpointer)
		for file in files:
			rtn.append(('graph', file))
	# return all file names
	return rtn


if __name__  == "__main__":
	print ("[Run]\tparameters detected: " + str(len(sys.argv)))
	if len(sys.argv) != 2:
		print ("[Run]\tusage: run <config_file>")
	else:
		cfg = config.getPipConfig(sys.argv[1])
		if cfg is not None:
			print('[Run]\tconfig loaded:')
			print(cfg)
			files = process(cfg)
		# write all files to 'output/log.txt'
		writeFiles(files)
