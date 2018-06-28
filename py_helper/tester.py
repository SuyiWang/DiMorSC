from visuals import visualcontainer
from visuals import volume
from visuals import graph
from dataio.ImageLoad import ImageLoader
from dataio.ImageWrite import ImageWriter
from PIL import Image
import numpy as np


# OP 4\6\7\8 is not consistent with others

IL = ImageLoader()

# print "##################"
# data1 = IL.LoadSingleTIF('Olfactory_Projection_Fibers/Image_Stacks/OP_1/1.tif')
# print data1.shape
# print "##################"

# data7 = IL.LoadSingleTIF('Olfactory_Projection_Fibers/Image_Stacks/OP_7/01.tif')
# print data7.shape

# print "##################"
# im = Image.open('Olfactory_Projection_Fibers/Image_Stacks/OP_7/01.tif')
# imarray7 = np.array(im)
# print imarray7.shape
# print np.amax(imarray7)

# print "op7 abs_difference " +str(np.amax(np.absolute(data7 - imarray7)))
#im.show()

data = IL.Load('Olfactory_Projection_Fibers/Image_Stacks/OP_7/', '.tif')

# for i in range(data.shape[-1]):
# 	IL.imagetestwrite(data[:,:,i], str(i));

# IL.imagetestwrite(data, '007.tif')

vc = visualcontainer()
vc.create_canvas()

volume.Draw(data, vc)
vc.Show()