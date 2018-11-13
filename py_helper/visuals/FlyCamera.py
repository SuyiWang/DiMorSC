from vispy import scene
from vispy.util.quaternion import Quaternion
import numpy as np
import math
from vispy.util import keys


# override several functions
class Fly(scene.cameras.FlyCamera):
    #   Modified from on_timer method
    #   Inertia in key board is removed
    #   center_loc = np.array(self._center, dtype='float32')
    #   self._center stores the center location


    # @property
    # def scale_factor(self):
    #   print 'scale_factor changed to' ,
    #   print self._scale_factor
    #   return self._scale_factor

    # here the _scale_factor is reset - hard coded
    # can be traced back to flycamera->perspectivecamera->baseCamera
    # it might be better if we set xlim first
    def __init__(self, **kwargs):
        super(Fly, self).__init__(**kwargs)
        self._timer.interval = 1.0/100

    def _set_range(self, init):
        """ Reset the view.
        """
        #PerspectiveCamera._set_range(self, init)
        # Stop moving
        self._speed *= 0.0

        # Get window size (and store factor now to sync with resizing)
        w, h = self._viewbox.size
        w, h = float(w), float(h)

        # Get range and translation for x and y
        x1, y1, z1 = self._xlim[0], self._ylim[0], self._zlim[0]
        x2, y2, z2 = self._xlim[1], self._ylim[1], self._zlim[1]
        rx, ry, rz = (x2 - x1), (y2 - y1), (z2 - z1)

        # Correct ranges for window size. Note that the window width
        # influences the x and y data range, while the height influences
        # the z data range.
        if w / h > 1:
            rx /= w / h
            ry /= w / h
        else:
            rz /= h / w

        # Do not convert to screen coordinates. This camera does not need
        # to fit everything on screen, but we need to estimate the scale
        # of the data in the scene.

        # Set scale, depending on data range. Initial speed is such that
        # the scene can be traversed in about three seconds.
        
        # this _scale_factor is not set appropriately
        #self._scale_factor = max(rx, ry, rz) / 3.0
        self._scale_factor = 500

        # Set initial position to a corner of the scene
        margin = np.mean([rx, ry, rz]) * 0.1
        self._center = x1 - margin, y1 - margin, z1 + margin
        self._center = (-100, -100, -50)

        # Determine initial view direction based on flip axis
        yaw = 45 * self._flip_factors[0]
        pitch = -90 - 20 * self._flip_factors[2]
        if self._flip_factors[1] < 0:
            yaw += 90 * np.sign(self._flip_factors[0])

        # Set orientation
        q1 = Quaternion.create_from_axis_angle(pitch*math.pi/180, 1, 0, 0)
        q2 = Quaternion.create_from_axis_angle(0*math.pi/180, 0, 1, 0)
        q3 = Quaternion.create_from_axis_angle(yaw*math.pi/180, 0, 0, 1)
        #
        self._rotation1 = (q1 * q2 * q3).normalize()
        self._rotation2 = Quaternion()

        # Update
        self.view_changed()

    def on_timer(self, event):
        """
        Timer event handler
        Parameters
        ----------
        event : instance of Event
            The event.
        """
        # Set relative speed and acceleration
        rel_speed = event.dt

        # Get what's forward
        pf, pr, pl, pu = self._get_directions()

        # Increase speed through acceleration
        # Rel acc removed!
        self._speed += self._acc

        # --- Determine new position from translation speed

        if self._speed[:3].any():

            # Create speed vectors, use scale_factor as a reference
            dv = np.array([1.0/d for d in self._flip_factors])
            #
            vf = pf * dv * rel_speed * self._scale_factor
            vr = pr * dv * rel_speed * self._scale_factor
            vu = pu * dv * rel_speed * self._scale_factor
            direction = vf, vr, vu

            # Set position
            center_loc = np.array(self._center, dtype='float32')
            center_loc += (self._speed[0] * direction[0] +
                           self._speed[1] * direction[1] +
                           self._speed[2] * direction[2])
            self._center = tuple(center_loc)

        # --- Determine new orientation from rotation speed

        roll_angle = 0

        # Calculate manual roll (from speed)
        if self._speed[3:].any():
            angleGain = np.array([1.0, 1.5, 1.0]) * 3 * math.pi / 180
            angles = self._speed[3:] * angleGain

            q1 = Quaternion.create_from_axis_angle(angles[0], -1, 0, 0)
            q2 = Quaternion.create_from_axis_angle(angles[1], 0, 1, 0)
            q3 = Quaternion.create_from_axis_angle(angles[2], 0, 0, -1)
            q = q1 * q2 * q3
            self._rotation1 = (q * self._rotation1).normalize()

        # Calculate auto-roll
        # NO AUTO-ROLL

        # Update
        if self._speed.any() or roll_angle or self._update_from_mouse:
            self._update_from_mouse = False
            self.view_changed()


        # Reduce speed. Simulate resistance. Using brakes slows down faster.
        # Note that the way that we reduce speed, allows for higher
        # speeds if keys ar bound to higher acc values (i.e. turbo)
        reduce = np.array([0.05, 0.05, 0.05, 0.1, 0.1, 0.1])
        reduce[self._brake > 0] = 0.2
        self._speed -= self._speed
        if np.abs(self._speed).max() < 0.05:
            self._speed *= 0.0

    def viewbox_mouse_event(self, event):
        """ViewBox mouse event handler
        Parameters
        ----------
        event : instance of Event
            The event.
        """
        if event.type == 'mouse_wheel':
            s = 1.1 ** - event.delta[1]
            self._scale_factor *= s

        if event.handled or not self.interactive:
            return

        if event.type == 'mouse_wheel':
            if not event.mouse_event.modifiers:
                # Move forward / backward
                #self._speed[0] += 0.5 * event.delta[1]
                pass
            elif keys.SHIFT in event.mouse_event.modifiers:
                # Speed
                s = 1.1 ** - event.delta[1]
                self.scale_factor /= s  # divide instead of multiply
                print('[FlyCamera] \tscale factor: %1.1f units/s' % self.scale_factor)
            return

        if event.type == 'mouse_press':
            event.handled = True

        if event.type == 'mouse_release':
            # Reset
            self._event_value = None
            # Apply rotation
            self._rotation1 = (self._rotation2 * self._rotation1).normalize()
            self._rotation2 = Quaternion()
            self._timer.stop()

        elif not self._timer.running:
            # Ensure the timer runs
            self._timer.start()

        if event.type == 'mouse_move':

            if event.press_event is None:
                return
            if not event.buttons:
                return

            # Prepare
            modifiers = event.mouse_event.modifiers
            pos1 = event.mouse_event.press_event.pos
            pos2 = event.mouse_event.pos
            w, h = self._viewbox.size

            if 1 in event.buttons and not modifiers:
                # rotate

                # get normalized delta values
                d_az = -float(pos2[0] - pos1[0]) / w
                d_el = +float(pos2[1] - pos1[1]) / h
                # Apply gain
                d_az *= - 0.5 * math.pi  # * self._speed_rot
                d_el *= + 0.5 * math.pi  # * self._speed_rot
                # Create temporary quaternions
                q_az = Quaternion.create_from_axis_angle(d_az, 0, 1, 0)
                q_el = Quaternion.create_from_axis_angle(d_el, 1, 0, 0)

                # Apply to global quaternion
                self._rotation2 = (q_el.normalize() * q_az).normalize()

            elif 2 in event.buttons and keys.CONTROL in modifiers:
                # zoom --> fov
                if self._event_value is None:
                    self._event_value = self._fov
                p1 = np.array(event.press_event.pos)[:2]
                p2 = np.array(event.pos)[:2]
                p1c = event.map_to_canvas(p1)[:2]
                p2c = event.map_to_canvas(p2)[:2]
                d = p2c - p1c
                fov = self._event_value * math.exp(-0.01*d[1])
                self._fov = min(90.0, max(10, fov))

        # Make transform be updated on the next timer tick.
        # By doing it at timer tick, we avoid shaky behavior
        self._update_from_mouse = True

    # camera controls
    def move(self, position):
        self._center = tuple(position)
        self.view_changed()
    
    '''
    // ref:https://www.khronos.org/opengl/wiki/GluLookAt_code
    float forward[3], side[3], up[3];
    float matrix2[16], resultMatrix[16];
    // --------------------
    forward[0] = center3D[0] - eyePosition3D[0];
    forward[1] = center3D[1] - eyePosition3D[1];
    forward[2] = center3D[2] - eyePosition3D[2];
    NormalizeVector(forward);
    // --------------------
    // Side = forward x up
    ComputeNormalOfPlane(side, forward, upVector3D);
    NormalizeVector(side);
    --------------------
    // Recompute up as: up = side x forward
    ComputeNormalOfPlane(up, side, forward);
    // --------------------
    matrix2[0] = side[0];
    matrix2[4] = side[1];
    matrix2[8] = side[2];
    matrix2[12] = 0.0;
    // --------------------
    matrix2[1] = up[0];
    matrix2[5] = up[1];
    matrix2[9] = up[2];
    matrix2[13] = 0.0;
    // --------------------
    matrix2[2] = -forward[0];
    matrix2[6] = -forward[1];
    matrix2[10] = -forward[2];
    matrix2[14] = 0.0;
    // --------------------
    matrix2[3] = matrix2[7] = matrix2[11] = 0.0;
    matrix2[15] = 1.0;
    // --------------------
    MultiplyMatrices4by4OpenGL_FLOAT(resultMatrix, matrix, matrix2);
    glhTranslatef2(resultMatrix,
                  -eyePosition3D[0], -eyePosition3D[1], -eyePosition3D[2]);
    // --------------------
    memcpy(matrix, resultMatrix, 16*sizeof(float));
    '''
    # one more example
    # https://github.com/vispy/vispy/issues/240
    def lookAt(self, viscenter):
        '''
        This version creates rotation quaterion that rotates front direction
        from global coordinates to eye coordinates.
        '''
        # compute actual rotate 
        front = np.array(viscenter) - np.array(self._center)
        front = normalize(front)
        
        front_in_cam = self.rotation.rotate_point(front)
        front_in_cam = normalize(front_in_cam)
        # in radius
        rollangle = math.acos(np.dot(front_in_cam, [0,0,-1]))
        rollaxis = np.cross(front_in_cam, [0,0,-1])
        rollaxis = normalize(rollaxis)

        quat = Quaternion.create_from_axis_angle(rollangle, 
                                                 rollaxis[0],
                                                 rollaxis[1],
                                                 rollaxis[2])
        # update rotation
        self.rotation = (quat * self.rotation).normalize()
        self.view_changed()


def normalize(v):
    norm = np.linalg.norm(v)
    if norm == 0: 
       return v
    return v / norm
