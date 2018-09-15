from PIL import Image
import numpy as np

'''
write 2d image file to png
'''
def write(A, filename, filetype = 'png'):
    #A = np.array(A, dtype=np.uint8)
    im = Image.fromarray(A)
    im.save(filename + '.' + filetype, filetype)