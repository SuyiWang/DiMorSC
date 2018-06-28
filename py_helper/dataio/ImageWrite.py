import numpy
import struct
from array import array
def ImageWriter(data, id, thd = -1, trans = [0, 0, 0]):
	newFile = open(id + "_dens.bin", "wb")
	# write to file
	double_array = array('d')
	counter = 0
	for (x,y,z),value in numpy.ndenumerate(data):
		if value > thd:
			counter+=1
			double_array.append(float(x) + trans[0])
			double_array.append(float(y) + trans[1])
			double_array.append(float(z) + trans[2])
			double_array.append(float(value))
	print "counted %d lines."%counter
	newFile.write(struct.pack('i', counter))
	double_array.tofile(newFile)
	newFile.close()
	