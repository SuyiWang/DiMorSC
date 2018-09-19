from itertools import cycle

from vispy import scene
from vispy.color import get_colormaps, BaseColormap
from vispy.visuals.transforms import STTransform
from vispy.scene.visuals import Text

from visuals.FlyCamera import Fly
from visuals import volume, graph, pointcloud

from dataio.struct import Data_Pointer
from dataio import imageLoader
from dataio import jsonHelper as json

import os

class MyCanvas():
    '''
    Add all visual object to this container
    '''

    def setup_canvas(self, canvas, targetview):
        # Create three cameras (Fly, Turntable and Arcball)
        fov = 60.
        cam1 = Fly(parent=targetview.scene, fov=fov, name='Fly')
        cam1._auto_roll = False
        cam2 = scene.cameras.TurntableCamera(parent=targetview.scene, fov=fov,
                                             name='Turntable', distance = 1000.0)
        cam3 = scene.cameras.ArcballCamera(parent=targetview.scene, fov=fov, name='Arcball', distance = 1000.0)
        targetview.camera = cam1  # Select Fly at first


        # Create an XYZAxis visual
        axis = scene.visuals.XYZAxis(parent=targetview)
        s = STTransform(translate=(50, 50), scale=(50, 50, 50, 1))
        affine = s.as_matrix()
        axis.transform = affine

        # create fps display
        self._fps = Text('', parent=self.view, color='white')
        self._fps.font_size = 8
        self._fps.pos = canvas.size[0]-100, 50

        # create colormaps that work well for translucent and additive volume rendering
        class TransFire(BaseColormap):
            glsl_map = """
            vec4 translucent_fire(float t) {
                return vec4(pow(t, 0.5), t, t*t, max(0, t*1.05 - 0.05));
            }
            """


        class TransGrays(BaseColormap):
            glsl_map = """
            vec4 translucent_grays(float t) {
                return vec4(t, t, t, t*0.05);
            }
            """

        # Setup colormap iterators
        opaque_cmaps = cycle(get_colormaps())
        translucent_cmaps = cycle([TransFire(), TransGrays()])
        opaque_cmap = next(opaque_cmaps)
        translucent_cmap = next(translucent_cmaps)


        # Implement axis connection with cam2
        @canvas.events.mouse_move.connect
        def on_mouse_move(event):
            if event.button == 1 and event.is_dragging:
                axis.transform.reset()

                axis.transform.rotate(cam2.roll, (0, 0, 1))
                axis.transform.rotate(cam2.elevation, (1, 0, 0))
                axis.transform.rotate(cam2.azimuth, (0, 1, 0))

                axis.transform.scale((50, 50, 0.001))
                axis.transform.translate((50., 50.))
                axis.update()


        # Implement key presses
        @canvas.events.key_press.connect
        def on_key_press(event):
            global opaque_cmap, translucent_cmap
            if event.text == '1':
                cam_toggle = {cam1: cam2, cam2: cam3, cam3: cam1}
                targetview.camera = cam_toggle.get(targetview.camera, cam2)
                print(targetview.camera.name + ' camera')
                if targetview.camera is cam2:
                    axis.visible = True
                else:
                    axis.visible = False
            elif event.text == '2':
                methods = ['mip', 'translucent', 'iso', 'additive']
                method = methods[(methods.index(volume1.method) + 1) % 4]
                print("Volume render method: %s" % method)
                cmap = opaque_cmap if method in ['mip', 'iso'] else translucent_cmap
                volume1.method = method
                volume1.cmap = cmap
            elif event.text == '4':
                if volume1.method in ['mip', 'iso']:
                    cmap = opaque_cmap = next(opaque_cmaps)
                else:
                    cmap = translucent_cmap = next(translucent_cmaps)
                volume1.cmap = cmap
            elif event.text == '0':
                cam1.set_range()
                cam3.set_range()
            elif event.text != '' and event.text in '[]':
                s = -0.025 if event.text == '[' else 0.025
                volume1.threshold += s
                th = volume1.threshold
                print("Isosurface threshold: %0.3f" % th)

        @canvas.events.resize.connect
        def resize(event=None):
            self._fps.pos = canvas.size[0]-100, 50


    def setup(self):
        # Prepare canvas
        self.canvas = scene.SceneCanvas(keys='interactive', size=(800, 600), show=True)
        self.canvas.measure_fps(callback = self.update_fps)

        # Set up a viewbox to display the image with interactive pan/zoom
        self.view = self.canvas.central_widget.add_view()
        
        self.setup_canvas(self.canvas, self.view)
        return


    def update_fps(self, toprint):
        printstr =  "%1.1f"%toprint
        self._fps.text = 'FPS: ' + printstr


    def insert(self, obj):
        # load all kinds of data, packaged as vobject
        ptr = load(obj.path)
        visual_item = _create_visual(ptr)
        if visual_item is not None:
            self._show(visual_item)

    def _show(self, vis):
        if not isinstance(vis, list):
            vis.parent = self.view.scene
        else:
            for item in vis:
                self._show(item)


def load(file, filetype = None):
    '''
    load data into Data_Pointer
    '''
    # known issue -> error if tries to load oldgraph 
    # and the path exists as a directory
    if os.path.isdir(file):
        return _loadFolder(file, filetype)
    else:
        return _loadFile(file, filetype)

def _loadFolder(foldername, filetype = None):
    '''
    returns a volume visual
    '''

    # Assuming foldername is a directory
    assert os.path.isdir(foldername), "[vobject]\tincorrect path detected"
    files = os.listdir(foldername)
    if len(files) == 0:
        print('[vobject]\tempty folder, nothing to do')
        return None
    
    # detect possible filetype
    # load majority file type (# >= 1)
    if filetype is None:
        filetype = imageLoader.majority_filetype(files)
    
    # load files, return None if load fail
    print("[vobject]\tFiletype to load: ", filetype)
    buff, downsample, orishape = imageLoader.loadVisFolder(foldername, filetype)
    data = {'volume': buff, 'type':'volume', 'shape':orishape}
    ptr = Data_Pointer(foldername, visdata=data, objtype='volume')
    
    return ptr


def _loadFile(path_id, filetype = None):
    '''
    returns a visual
    '''

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
    
    ptr = _getdata(path_id, filetype)
    return ptr

def _getdata(path_id, filetype):
    '''
    create Data_Pointer containing loaded visdata
    '''

    # use different function according to filetype
    if filetype == 'tif':   
        visdata, downsample, orishape = imageLoader.loadVisVolume(path_id, filetype)
        data = {'volume':visdata, 'type':'volume'}
        ptr = Data_Pointer(path_id, visdata=data, objtype='volume')
        return ptr
    elif filetype == 'oldgraph':
        edge, vert = graph.Load(path_id)
        if edge is not None and vert is not None:
            data = {'vert':vert, 'edge':edge, 'type':'graph'}
            ptr = Data_Pointer(path_id, visdata=data, objtype='graph')
            return ptr
        else:
            print('[vobject]:\tskipped empty file', path_id)
            return None
    elif filetype == 'sc':
        # simplicial complex is visualized as graph
        # and the triangle information is ignored
        edge, vert = graph.LoadSimplexGraph(path_id)
        if edge is not None and vert is not None:
            data = {'vert':vert, 'edge':edge, 'type':'graph'}
            ptr = Data_Pointer(path_id, visdata=data, objtype='graph')
            return ptr
        else:
            print('[vobject]:\tskipped empty file', path_id)
            return None
    elif filetype == 'ini':
        edge, vert = graph.loadini(path_id)
        if edge is not None and vert is not None:
            data = {'vert':vert, 'edge':edge, 'type':'graph'}
            ptr = Data_Pointer(path_id, visdata=data, objtype='graph')
            return ptr
        else:
            print('[vobject]:\tskipped empty file', path_id)
            return None
    elif filetype == 'swc':
        edge, vert = graph.loadSWC(path_id)
        if edge is not None and vert is not None:
            data = {'vert':vert, 'edge':edge, 'type':'graph'}
            ptr = Data_Pointer(path_id, visdata=data, objtype='graph')
            return ptr
        else:
            print('[vobject]:\tskipped empty file', path_id)
            return None
    elif filetype == 'dens':
        vert, downsample = pointcloud.load(path_id)
        data = {'vert':vert, 'type':'cloud'}
        return Data_Pointer(path_id, visdata=data, objtype='cloud')
    elif filetype == 'json':
        # get all files
        inputlist = json.parse_data(path_id, jsontype='block')
        findpath = os.path.dirname(path_id) + '/'
        # load batch
        ptrlist = []
        for obj in inputlist:
            loaded = load(findpath+obj['file'])
            if loaded is not None:
                loaded.offset = obj['pos'][0:3]
                #loaded.offset = obj['pos']
                ptrlist.append(loaded)
        
        if len(ptrlist) != 0:
            # force visdata contains a list of vobjects
            blockobj = Data_Pointer(path_id, 
                               objtype='json',
                               visdata={'vollist':ptrlist})
            return blockobj
        else:
            return None
    else:
        return None

def _create_visual(ptr):
    # create vispy visual using _visdata
    
    if ptr is None:
        return None

    if ptr.type == 'json':
        return volume.multi_draw(ptr._visdata['vollist'])
    elif ptr.type == 'cloud':
        return pointcloud.draw(ptr._visdata['vert'])
    elif ptr.type == 'graph':
        return graph.draw(ptr._visdata['edge'],
                          ptr._visdata['vert'])
    elif ptr.type == 'volume':
        return volume.draw(ptr._visdata['volume'])
    else:
        return None
