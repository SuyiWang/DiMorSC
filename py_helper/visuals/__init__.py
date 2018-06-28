from itertools import cycle

import numpy as np

from vispy import app
from vispy import scene

from vispy.color import get_colormaps, BaseColormap
from vispy.visuals.transforms import STTransform
from FlyCamera import Fly

class visualcontainer():
    '''
    Add all visual object to this container
    Use clear to clear drawing
    
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

    def create_canvas(self):
    	# Prepare canvas
        self.canvas = scene.SceneCanvas(keys='interactive', size=(800, 600), show=True)
        self.canvas.measure_fps()


        # Set up a viewbox to display the image with interactive pan/zoom
        self.view = self.canvas.central_widget.add_view()

        self.visuallist = []
        
        self.setup_canvas(self.canvas, self.view)
        return
    def clear(self):
        for v in self.visuallist:
            v.parent = None
            del v
        self.visuallist = []

    def Show(self):
        self.canvas.show()
        app.run()