from vispy import scene
import numpy as np
import math


#	Todo: switch to mouse keyboard driven event. not timer
class Fly(scene.cameras.FlyCamera):
	#   Override this function
	#	Inertia in key board is removed
	#   In the future, this can be controlled by pure keyboard event instead 
	#   of time event

	#   center_loc = np.array(self._center, dtype='float32')
	#   self._center stores the center location
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