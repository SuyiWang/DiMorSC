
from dataio import splitter
import unittest

import warnings
def ignore_warnings(test_func):
    def do_test(self, *args, **kwargs):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            test_func(self, *args, **kwargs)
    return do_test

# python3 -m unittest discover py_helper/
class test_upper_bound(unittest.TestCase):
	'''
	test binary search algorithm
	'''
	def setUp(self):
		self.data = [0, 0, 1, 1, 2, 3, 4, 5]
		self.func = lambda x: self.data[x]
		self.func.search_range = (0, len(self.data))

	@ignore_warnings
	def runTest(self):
		assert splitter._upper_bound(self.func, 0) == 1
		assert splitter._upper_bound(self.func, 1) == 3
		assert splitter._upper_bound(self.func, 5) == 7
		assert splitter._upper_bound(self.func, 10) == 7
		assert splitter._upper_bound(self.func, -1) == -1
		assert splitter._upper_bound(self.func, 1.5) == 3

	def tearDown(self):
		pass

class test_lower_bound(unittest.TestCase):
	'''
	test binary search algorithm
	'''
	def setUp(self):
		self.data = [0, 0, 1, 1, 2, 3, 4, 5]
		self.func = lambda x: self.data[x]
		self.func.search_range = (0, len(self.data))

	@ignore_warnings
	def runTest(self):
		assert splitter._lower_bound(self.func, 0) == 0
		assert splitter._lower_bound(self.func, 1) == 2
		assert splitter._lower_bound(self.func, 5) == 7
		assert splitter._lower_bound(self.func, 10) == 8
		assert splitter._lower_bound(self.func, -1) == 0
		assert splitter._lower_bound(self.func, 1.5) == 4

	def tearDown(self):
		pass