'''
Module for loading all types of image files

Methods:
    (numpy.ndarray) loadFolder(string path, string filetype):
        Loads volume with given extension from a folder
        Try to load all file if filetype is not specified

    (numpy.ndarray, bool) loadVisFolder(string path, string filetype):
        Loads volume with given extension from a folder
        And DOWNSAMPLE it to visualization size
        return True if downsampled

    (numpy.ndarray) load(string path, string suffix):
        Loads data from a file or folder

    (list[string]) getfilelist(string path, string suffix):
        acquire file list in numeric order, if there is number in filename

    (numpy.ndarray) loadVolume(string filename, string filetype = None):
        Reads a volume file.

    (numpy.ndarray, bool) loadVisVolume(string filename, string filetype = None):
        Reads a volume file and DOWNSAMPLE it.
        return True if downsampled

    (numpy.ndarray) loadImage(string filename, string filetype = None):
        Reads a 2D image.

    (None) imageTestWrite(numpy.ndarray A, string filename)
        Writes a 2d image to file.

    [internal methods]
    tryint:
    alphanum_key:
        function for sorting files according to split string and number
'''
#[TODO] downsample factor now hard coded.


import numpy as np
import scipy.ndimage as ndimage
from sys import argv
import glob, os
# slow import - but fine for now
from skimage import io
from skimage.measure import block_reduce
import re


# sometimes need this, esp when the following error message is shown
# PIL.Image.DecompressionBombError: Image size (352716000 pixels) 
# exceeds limit of 178956970 pixels, could be decompression bomb DOS attack.
from PIL import Image
Image.MAX_IMAGE_PIXELS = 352716000


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


'''
loads image files from a folder
filetype decides which subset to load
filetype == '' means all files

(string) path, (string) filename
'''
def loadFolder(path, filetype = None):
    if path[-1] is not '/':
        path = path + '/'
    # get all files of "filetype"
    files = GetFilelist(path, filetype)
    data = []
    for file in files:
        # filetype tells which module to use for loading data
        data.append(loadImage(path+file))

    # np.stack is fastest for axis=0
    # locality efficient. so we stick to axis=0
    # z-axis is the first coordinate, then x,y
    # this should made consistent through whole program
    rtn = np.stack(data)
    # skip an empty line
    print("    [ImageLoader]\tFiles stacked with size", rtn.shape)
    return rtn


'''
loadFolder method with downsample.
'''
def loadVisFolder(path, filetype = None):
    if path[-1] is not '/':
        path = path + '/'
    # get all files of "filetype"
    files = GetFilelist(path, filetype)
    print("    [ImageLoader]\tTotal %d files detected"%len(files))
    data = []
    orishape = None

    # Test size
    # If size > 1024*1024, down sample
    # if length > 500, down sample
    buff = loadImage(path+files[0])
    limit = []
    downsample = False
    for dima, dimb in zip(buff.shape, [512, 512]):
        limit.append((dima-1) // dimb + 1)
        if (limit[-1] > 1): downsample = True
    limit = tuple(limit)
    orishape = (len(files), ) + buff.shape

    for file in files:
        # filetype tells which module to use for loading data
        buff = loadImage(path+file)
        if downsample:
            buff = block_reduce(buff, block_size=limit, func=np.mean)
        data.append(buff)

    # np.stack is fastest for axis=0
    # locality efficient. so we stick to axis=0
    # z-axis is the first coordinate, then x,y
    data = np.stack(data)

    # accordingly, this down sample is done on limit[0] instead of 
    # limit[-1] as originally written
    limit = ((len(files) - 1)//512 + 1, ) + limit
    if (limit[0] > 1): 
        downsample = True
        data = block_reduce(data, block_size=(limit[0], 1, 1), func = np.mean)

    # skip an empty line
    if downsample:
        print("    [ImageLoader]\tDownsample factor", limit)
    print("    [ImageLoader]\tFiles stacked with size", data.shape)
    # make sure correct data is loaded later
    return data, downsample, orishape


'''
Loads all types of images: Single and folders
if folder, suffix is 'tif' unless specified.
'''
def load(filepath, suffix=None):
    if os.path.isdir(filepath):
        return loadFolder(filepath, suffix)
    elif os.path.isfile(filepath):
        #  This should be changed to consider more file types
        data = loadVolume(filepath)
        return data
    else:
        print("    [ImageLoader]\tUnable to load %s" %filepath)
        return None


def majority_filetype(files):
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
    return filetype


# TODO: select only files without .ext if suffix == ''
'''
get file list in numeric order
(string) path, (string) suffix
'''
def GetFilelist(path, suffix):
    print("    [ImageLoader]\tReading files in", path)
    selected = []
    files = os.listdir(path)
    if suffix is None:
        suffix = majority_filetype(files)
    pat = re.compile('.*\.'+suffix+'$', re.I)
    for file in files:
        if pat.match(file) is not None:
            selected.append(file)
    return sorted(selected, key=alphanum_key)


'''
Load 3d Volume file, e.g. tiff

(string) filename, (string) filetype
'''
def loadVolume(filename, filetype = None):
    if filename is None:
        print ("    [ImageLoader]\tNothing to do.")
        return None
    printfilename = filename.split('/')[-1]
    im  = io.imread(filename)
    if im is None:
        print ("    [ImageLoader]\tFailed to load %s"%filename)
        return None
    #  [TODO?] coule be changed for better precision
    #  F8 is 64bit = double

    # if im.dtype != np.uint8:
    #   maxx = np.amax(im)
    #   #print "max value: " + str(maxx)
    #   rtn = np.array(im, dtype='f8')
    #   rtn = rtn*(255.0/maxx)
    #   ###Vispy has no constraints for uint8 or uint16

    # else:
    #   rtn = np.array(im, dtype='f8')
    rtn = im
    
    msg = "    [ImageLoader]\tImage %s loaded. Size: "%printfilename  + str(rtn.shape)
    print(msg)
    return rtn


'''
Load 3d Volume file and downsample, e.g. tiff

(string) filename, (string) filetype
'''
def loadVisVolume(filename, filetype = None):
    orishape = None
    if filename is None:
        print ("[ImageLoader]\tNothing to do.")
        return None, False, None
    printfilename = filename.split('/')[-1]
    im  = io.imread(filename)
    if im is None:
        print ("[ImageLoader]\tFailed to load %s"%filename)
        return None, False, None
    
    msg = "    [ImageLoader]\tImage %s loaded. Size: "%printfilename  + str(im.shape)
    print(msg)
    orishape = im.shape

    # If size > 500*1024*1024 reshape
    limit = []
    downsample = False
    for dima, dimb in zip(im.shape, [500, 512, 512]):
        limit.append((dima-1) // dimb + 1)
        if (limit[-1] > 1): downsample = True
    limit = tuple(limit)

    if downsample:
        im = block_reduce(im, block_size=tuple(limit), func=np.mean)
        print('    reshaped to ', im.shape)
    # must make sure correct data is loaded later
    return im, downsample, orishape


'''
Load a 2D image, will compress a third dimension with sum()
if the file has one.

filetype now is not used yet

(string) filename, (string) filetype
'''
def loadImage(filename, filetype = None):
    if filename is None:
        print ("    [ImageLoader]\tNothing to do.")
        return None
    printfilename = filename.split('/')[-1]
    im  = io.imread(filename)
    if im is None:
        print ("    [ImageLoader]\tFailed to load %s"%filename)
        return None
    
    if (im.ndim > 2):
        im = np.sum(im, axis=2)
    # rescale grayscale value might cause several problems
    # removed for now
    msg = "    [ImageLoader] Image %s loaded. Size: "%printfilename  + str(im.shape)
    print(msg, "Max Pixel:", (np.amax(im)))
    return im


'''
write 2d image file to png
'''
def imageTestWrite(A, filename, filetype = 'PNG'):
    A = np.array(A, dtype=np.uint8)
    im = Image.fromarray(A)
    file = filename.split('/')[-1]
    file = file.split('.')[0]
    im.save('test/'+file+'.png','PNG')


if __name__  == "__main__":
    print ("[ImageLoader]\tparameters detected: " + str(len(argv)))
    if argv[2] == '--tif':
        data = load(argv[1], 'tif')
    else:
        print ("[ImageLoader]\tusage: ImageLoad <filename> <filetype>")
