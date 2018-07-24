from os import path


template_pip = """# Comments
# Do not remove unused parameters
# [Required for all] input filename
{inputname}
# [Required for all] file type volume/cloud/sc/graph
{type}
# [Required for volume] Smoothing Sigma.
{sigma}
# [Required for volume] Weak threshold.
{densThd}
# DiMorSC parameters.
# [Required for DiMorSC] Persistence threshold
{persist}
# [Required for graph2tree] Root
{root}
# [Required for graph2tree] Saddle threshold
{saddleThd}

"""

def getPipConfig(file):
	rtn = []
	with open(file, "r") as fp:
		line = fp.readline()
		while line:
			if line.startswith('#'):
				line = fp.readline()
				continue
			# readline has \n in the end
			rtn.append(line.rstrip('\n'))
			line = fp.readline()
	fp.close()

	if len(rtn) != 7:
		print('[getPipConfig] Exactly 6 parameters must be provided. doing nothing.')
		return None
	else:
		if not path.isdir(rtn[0]) and not path.isfile(rtn[0]):
			print('[getPipConfig]\t', rtn[0], 'does not exist')
			return None
		
		# TODO
		# check sigma should be positive
		# check filetype should be volume/cloud/sc/graph
		# check densThd should be number
		# check persist should be positive
		# check root should be a list of length 3
		# check saddle Thd should be number
		return rtn

def getDefaultPip():
	config = [
	'data/Olfactory_OP7_trunc.tif',
	'volume',
	'2',
	'0.1',
	'5',
	'0 0 0',
	'0'
	]
	return config

def writePipConfig(filename, pip_config):
	f = open(filename, 'w')
	if f is not None:
		to_write = template_pip.format(
			inputname = pip_config[0],
			type = pip_config[1],
			sigma = pip_config[2],
			densThd = pip_config[3],
			persist = pip_config[4],
			root = pip_config[5],
			saddleThd = pip_config[6],
			)
		f.write(to_write)
		f.close()
	else:
		print('[config]\tError creating file:', filename)


template_graph = """# vertex filename
{vertname}
# edge filename
{edgename}
# output prefix: program writes to output/0.swc
{outname}
# root
{root}
# threshold for simplification
{sadd_thd}
"""

def getGraphConfig(file):
	print('[config]\tget graph config: TODO')

def writeGraphConfig(filename, graph_config):
	f = open(filename, 'w')
	if f is not None:
		to_write = template_graph.format(
			vertname=graph_config[0],
			edgename=graph_config[1],
			outname=graph_config[2],
			root=graph_config[3],
			sadd_thd=graph_config[4],
			)
		f.write(to_write)
		f.close()
	else:
		print('[config]\tError creating file:', filename)

def getDefaultGraph():
	config = [
	'0_vert.txt',
	'0_edge.txt',
	'output/0',
	'0 0 0',
	'0'
	]
	return config

if __name__  == "__main__":
	config = getDefaultGraph()
	writeGraphConfig('output/test.ini', config)
	
