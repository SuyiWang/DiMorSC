# Split image stacks to smaller cubes

# a call looks like this:
# file_info, offset = splitter.to_blocks(ptr, size)
# ptr: Data_Pointer
# size: string, specification of block size
import os
from dataio import imageLoader
from dataio import imageWriter
from math import ceil
from dataio import jsonHelper as json


def to_blocks(ptr, size, overlap = [3, 3, 3], workpath='output/'):
    if not os.path.isdir(ptr._path):
        return split_3D_block(ptr, size, overlap, workpath=workpath)
    else:
        return split_picwise_block(ptr, size, overlap, workpath=workpath)


def split_3D_block(ptr, size, overlap, workpath='output/'):
    # load data
    ptr._visdata['volume'] = imageLoader.load(ptr.path)
    ptr._visdata['shape'] = ptr._visdata['volume'].shape
    # shape in ZYX order
    file_struct = cal_split(ptr._visdata['shape'], size, overlap)

    # split ndarray to desired blocks
    # get folder name
    # write for each block
    data_name = ptr._get_name()
    base_folder = workpath + data_name + '_split/'
    #if folder does not exist, create
    os.makedirs(base_folder, exist_ok=True)
    
    wholedata = ptr._visdata['volume']
    for filenum in range(len(file_struct)):
        write_3D_block(wholedata, 
                       base_folder + str(file_struct[filenum]['file']),
                       file_struct[filenum]['pos'])
    

    # finishing touch, write json describing the structure of blocks
    jsondata = {"type":"block", "data":file_struct, "overlap":overlap}
    jsonName = base_folder + data_name +'.json'
    json.write(jsondata, jsonName)
    return jsonName


def write_3D_block(data, foldername, offset):
    # if folder does not exist, create it
    os.makedirs(foldername, exist_ok=True)
    if foldername[-1] != '/':
        foldername = foldername + '/'

    layernum = 0
    # offset should be guaranteed to be 3D
    # loop over all layers, data[z, y, x]. z corresponds to vertical layers
    for i in range(offset[0], offset[3]):
        imageWriter.write(data[i, offset[1]:offset[4], offset[2]:offset[5]], 
                          foldername + str(layernum),
                          'png')
        layernum += 1

def split_picwise_block(ptr, size, overlap, workpath='output/'):
    # generate folder name
    data_name = ptr._get_name()
    data_folder = workpath+data_name+'_split/'
    os.makedirs(data_folder, exist_ok=True)

    # get file number from path
    files = imageLoader.GetFilelist(ptr.path)
    volume_shape = [len(files), ]

    # z loop
    # for each image
    zaxis = 0
    for file in files:
        # Assuming no missing files

        # filetype tells which module to use for loading data
        img = imageLoader.loadImage(ptr.path+'/'+file)

        if len(volume_shape) != 3:
            # shape[0] should be y axis, shape[1] should be x.
            volume_shape.append(img.shape[0])
            volume_shape.append(img.shape[1])
            file_struct = cal_split(volume_shape, size, overlap)

        # file_struct must exist now.
        # find upper/lower bound for z axis.
        bounds = _find_bounds(file_struct, zaxis)

        # get folder name
        # loop over all structures, a folder is created for each struct
        # filestruct in this loop is guaranteed to intersect with give pic
        for foldernum in range(bounds[0], bounds[1]):
            struct = file_struct[foldernum]
            offset = struct['pos']
            # write one file
            _write_2D_img(img[offset[1]:offset[4], offset[2]:offset[5]], 
                          data_folder + struct['file'],
                          zaxis)

        zaxis += 1

    # finishing touch, write json describing the structure of blocks
    jsondata = {"type":"block", "data":file_struct, "overlap":overlap}
    jsonName = data_folder + data_name +'.json'
    json.write(jsondata, jsonName)
    return jsonName

def _find_bounds(file_struct, z):
    # it's binary search
    # struct:file_struct is offset zl = struct[0] zu = struct[3]
    
    # lower bound L is the last x: z_upper(x) <= z
    # upper bound U is the first x; z_lower(x) >= z
    z_upper = lambda x: file_struct[x]['pos'][3]
    z_upper.search_range = (0, len(file_struct))

    z_lower = lambda x: file_struct[x]['pos'][0]
    z_lower.search_range = (0, len(file_struct))
    
    L = _lower_bound(z_upper, z)
    U = _upper_bound(z_lower, z)

    # solution always exist
    return (L, U+1)

def _upper_bound(array_func, target):
    '''
    suppose f(x) monotonic inc
    returns the last index x where f(x) <= target
    '''
    # array_obj should have attribute search_range()
    # array_obj(x)

    # [left, right)
    f = array_func
    left, right = f.search_range

    while left<right:
        mid = (left + right) // 2
        if f(mid) > target:
            right = mid
        else:
            left = mid + 1
    return right - 1

def _lower_bound(array_func, target):
    '''
    suppose f(x) monotonic inc
    returns the first index x where f(x) >= target
    '''

    f = array_func
    left, right = f.search_range

    while left < right:
        mid = (left + right) // 2
        if f(mid) >= target:
            right = mid
        else:
            left = mid + 1
    return right


def _write_2D_img(imgpiece, folder, filenum):
    os.makedirs(folder, exist_ok=True)
    if folder[-1] != '/':
        folder = folder + '/'
    imageWriter.write(imgpiece, folder + str(filenum), 'png')

def cal_split(datashape, targetsize, overlap):
    '''
    calculate actual shape for each block
    output:
        a list of dictionaries, each containing key words "file" and "pos"
        The output is JSON ready
    '''

    blksize = targetsize.split(' ')
    DIM = 3
    assert len(datashape)==DIM or len(blksize)==DIM or len(overlap)==DIM, \
           "volume, target block or overlap dimension mismatch"

    # calculate cuts
    # for every dimension: L is total pixel, a is target size, x is overlap
    # output: (i, [i*(a-x), i*(a-x) + a-1]), where i in [0, ceil((L-a)/(a-x))]
    # [] means inclusive
    intsize = []
    bounds = []
    for i in range(DIM):
        if blksize[i] == '*':
            intsize.append(datashape[i])
            bounds.append(1)
        else:
            intsize.append(int(blksize[i]))
            bounds.append( 
                ceil( (datashape[i]-intsize[i]) / (intsize[i]-overlap[i]) ) + 1
                )
    
    # This iteration is in ZYX order.
    # so generated trans/'pos'/offset are in ZYX order
    # where YX is a plane of image.
    file_struct = []
    for i in range(bounds[0]):
        for j in range(bounds[1]):
            for k in range(bounds[2]):
                namestr = str(serialize(i, j, k, bounds))
                ir = get_range(i, intsize[0], overlap[0], datashape[0])
                jr = get_range(j, intsize[1], overlap[1], datashape[1])
                kr = get_range(k, intsize[2], overlap[2], datashape[2])
                offset = (ir[0], jr[0], kr[0], ir[1], jr[1], kr[1])
                file_struct.append({'file':namestr, 'pos':offset})
    # then for 3 dimension we have LWH for maximum number blocks.
    # convert each i,j,k to a canonical index c in range 0-L, 0-W, 0-H
    # c is the final filenumber. c -> i,j,k -> offset(i,j,k), given LWH
    return file_struct

def serialize(i, j, k, bounds):
    return i*bounds[1]*bounds[2] + j*bounds[2] + k

def deserialize(idx, bounds):
    k = idx % bounds[2]
    idx = idx/bounds[2]
    j = idx % bounds[1]
    idx = idx/bounds[1]
    i = idx
    return i, j, k

def get_range(i, a, x, limit):
    # return value [), compatible with range function
    return (i*(a-x), min(i*(a-x) + a, limit))