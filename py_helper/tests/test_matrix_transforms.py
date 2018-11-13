from visuals import FlyCamera
import numpy as np
import unittest
from vispy.visuals.transforms import MatrixTransform
import math
from pyquaternion import Quaternion

import warnings
def ignore_warnings(test_func):
    def do_test(self, *args, **kwargs):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            test_func(self, *args, **kwargs)
    return do_test

# python3 -m unittest discover py_helper/
class test_transformation(unittest.TestCase):
	'''
	test binary search algorithm
	'''
	def setUp(self):
		# points are given in origin/front/right/up order
		# pglobal are mapped to peye
		self.peye = np.array([[0,0,0], [0,0,-1], [1,0,0], [0,1,0]])
		self.pglobal = np.array([[2,2,2], [2,3,2], [2,2,3], [1,2,2]])
		self.protate = np.array([[0,0,0], [2,3,2], [2,2,3], [1,2,2]])

	@ignore_warnings
	def runTest(self):
		# test affine mapping
		### Case 1
		cam = FlyCamera.Fly()
		M = cam.transform
		# M maps pglobal to peye
		M.set_mapping(self.pglobal, self.peye)
		self.assertTrue(np.allclose(self.peye, M.map(self.pglobal)[:,:3]))

		### Case 2
		M.set_mapping(self.protate, self.peye)
		self.assertTrue(np.allclose(self.peye, M.map(self.protate)[:,:3]))
		

	def tearDown(self):
		pass