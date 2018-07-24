'''
[CLASS] vobject
a data portal used in front end 
and _path stores real data pointer for backend.

Properties:
	(numpy.ndarray) _visdata:
		contains data for visualization. Memory consuming
	(string) _type:
		"graph/sc/volume/cloud"
		sc is visualized as graph
	(string) _path:
		absolute path of data source
	(int) _ID:
		id for passing information
	(bool) _dataloaded:
		True if the data is completely loaded
	(string) _modtime:
		Last modification time and date of the loaded data
Methods:
	(string) get_default_nickname():
		get name.ext from path
	(string) get_name():
		get name without extension
	(void) setID():
	(bool) update():
		refetch the data for updates
'''


from dataio import imageLoader
import os
from visuals import graph, pointcloud
import glob


class vobject():
	# [Types] string:path, numpy.ndarray:visdata, string:objtype, 
	#		  int:ID, string:modtime
	# objtype: graph/cloud/volume/sc
	def __init__(self, path, visdata=None, objtype=None,
				 ID=None, modtime = None, complete = True):
		# path should be FULL PATH/abspath
		self._path = path
		self._visdata = visdata
		self._type = objtype
		self._ID = ID
		self._modtime = modtime
		# if dataloaded, use _visdata
		# otherwise, load from _path
		self._dataloaded = complete

	def get_default_nickname(self):
		# in windows does this cause problem?
		return self._path.split('/')[-1]

	def get_name(self):
		fullname = self.get_default_nickname()
		namesplit = fullname.split('.')
		if len(namesplit) > 1:
			return namesplit[-2]
		else:
			return namesplit[-1]

	def setID(self, ID):
		self._ID = ID

	def update(self):
		newobj = load(self._path)
		self._visdata = newobj._visdata
		self._modtime = newobj._modtime
		self._dataloaded = newobj._dataloaded
		# check type and dataloaded
		del newobj


# load all kinds of data, packaged as vobject
def load(file, filetype = None):
	if os.path.isdir(file):
		return _loadFolder(file, filetype)
	else:
		return _loadFile(file, filetype)


def _loadFolder(foldername, filetype = None):
	# Assuming foldername is a directory
	assert os.path.isdir(foldername), "[vobject]\tincorrect path detected"
	files = os.path.listdir(foldername)
	if len(files) == 0:
		print('[vobject]\tempty folder, nothing to do')
		return None
	
	# detect possible filetype
	# load majority file type (# >= 1)
	if filetype is None:
		suffices = {}
		maxcount = 0
		for file in files:
			splitted = file.split('.')
			if len(splitted) > 1:
				ext = splitted[-1].lower()
			else:
				ext = ''

			if ext in suffices:
				cnt = suffices[ext] + 1
				suffices[ext] = cnt
			else:
				cnt = 1
				suffices[ext] = 1
			# print(ext, suffices[ext])
			if cnt > maxcount:
				maxcount = cnt

		for key, value in suffices.items():
			if value == maxcount:
				filetype = key
				break
	
	# load files, return None if load fail
	print("[vobject]\tFiletype to load: ", filetype)
	buff, downsample = imageLoader.loadVisFolder(foldername, filetype)
	data = {'volume': buff, 'type':'volume'}
	vobj = vobject(foldername, data, 'volume', complete=not downsample)
	return vobj


def _loadFile(path_id, filetype = None):
	# detect file type if not specified
	if filetype is None:
		splitfile = path_id.split('.') 
		if len(splitfile) == 1:
			# no suffix is graph
			filetype = 'oldgraph'
			print ('[vobject]\tloading old graph format: ', path_id)
		else:
			filetype = splitfile[-1]
			print ('[vobject]\tloading ', filetype, ': ', path_id)
			assert os.path.isfile(path_id), "[vobject]\tFile does not exist"
	
	# use different function according to filetype
	if filetype == 'tif':	
		visdata, downsample = imageLoader.loadVisVolume(path_id, filetype)
		data = {'volume':visdata, 'type':'volume'}
		vobj = vobject(path_id, data, 'volume', complete=not downsample)
		return vobj
	elif filetype == 'oldgraph':
		edge, vert = graph.Load(path_id)
		if edge is not None and vert is not None:
			data = {'vert':vert, 'edge':edge, 'type':'graph'}
			vobj = vobject(path_id, data, 'graph')
			return vobj
		else:
			print('[vobject]:\tskipped empty file', path_id)
			return None
	elif filetype == 'sc':
		# simplicial complex is visualized as graph
		# and the triangle information is ignored
		edge, vert = graph.LoadSimplexGraph(path_id)
		if edge is not None and vert is not None:
			data = {'vert':vert, 'edge':edge, 'type':'graph'}
			return vobject(path_id, data, 'sc')
		else:
			print('[vobject]:\tskipped empty file', path_id)
			return None
	elif filetype == 'ini':
		edge, vert = graph.loadini(path_id)
		if edge is not None and vert is not None:
			data = {'vert':vert, 'edge':edge, 'type':'graph'}
			return vobject(path_id, data, 'graph')
		else:
			print('[vobject]:\tskipped empty file', path_id)
			return None
	elif filetype == 'swc':
		edge, vert = graph.loadSWC(path_id)
		if edge is not None and vert is not None:
			data = {'vert':vert, 'edge':edge, 'type':'graph'}
			return vobject(path_id, data, 'graph')
		else:
			print('[vobject]:\tskipped empty file', path_id)
			return None
	elif filetype == 'dens':
		vert, downsample = pointcloud.load(path_id)
		data = {'vert':vert, 'type':'cloud'}
		return vobject(path_id, data, 'cloud', complete=not downsample)
	else:
		return None
