from vispy import scene
from vispy.util.quaternion import Quaternion
import numpy as np
import math
from vispy.util import keys

# override several functions
class Fly(scene.cameras.FlyCamera):
	#   Modified from on_timer method
	#	Inertia in key board is removed
	#   center_loc = np.array(self._center, dtype='float32')
	#   self._center stores the center location


	# @property
	# def scale_factor(self):
	# 	print 'scale_factor changed to' ,
	# 	print self._scale_factor
	# 	return self._scale_factor

	# here the _scale_factor is reset - hard coded
	# can be traced back to flycamera->perspectivecamera->baseCamera
	# it might be better if we set xlim first
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
		#print 'updating'
		rel_speed = event.dt
		rel_acc = 0.5

		# Get what's forward
		pf, pr, pl, pu = self._get_directions()

		# Increase speed through acceleration
		# Note that self._speed is relative. We can balance rel_acc and
		# rel_speed to get a nice smooth or direct control
		self._speed += self._acc * rel_acc

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
		if self.auto_roll:
			up = {'x': (1, 0, 0), 'y': (0, 1, 0), 'z': (0, 0, 1)}[self.up[1]]
			up = np.array(up) * {'+': +1, '-': -1}[self.up[0]]

			def angle(p1, p2):
				return np.arccos(p1.dot(p2))
			#au = angle(pu, (0, 0, 1))
			ar = angle(pr, up)
			al = angle(pl, up)
			af = angle(pf, up)
			# Roll angle that's off from being leveled (in unit strength)
			roll_angle = math.sin(0.5*(al - ar))
			# Correct for pitch
			roll_angle *= abs(math.sin(af))  # abs(math.sin(au))
			if abs(roll_angle) < 0.05:
				roll_angle = 0
			if roll_angle:
				# Correct to soften the force at 90 degree angle
				roll_angle = np.sign(roll_angle) * np.abs(roll_angle)**0.5
				# Get correction for this iteration and apply
				angle_correction = 1.0 * roll_angle * math.pi / 180
				q = Quaternion.create_from_axis_angle(angle_correction,
													  0, 0, 1)
				self._rotation1 = (q * self._rotation1).normalize()

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