#!/usr/bin/python3
import unittest
import argparse
from pathlib import Path
print('Running' if __name__ == '__main__' else 'Importing', Path(__file__).resolve())

#sys.path.append('/home/oliver/PycharmProjects/Project_scripts/refactored_code/existing_code')
from onehelixres_refactored import single_helix_parser

class TestSum(unittest.TestCase):

    def test_sum(self):
        self.assertEqual(sum([1, 2, 3]), 6, "Should be 6")

    def test_sum_tuple(self):
        self.assertEqual(sum((1, 2, 2)), 6, "Should be 6")

class TestSum(unittest.TestCase):
    def test_list_int(self):
        """
        Test that it can sum a list of integers
        """
        data = [1, 2, 3]
        result = sum(data)
        self.assertEqual(result, 6)

# class Testonehelixres(unittest.TestCase):
#     def test_onehliex(self):
#         """
#         Tests that onehelixres.py output the right bloody thing
#         """
#         with open('1ct5_old.1hr', 'r') as old_file:
#             with open('1ct5.1hr', 'r') as new_file:
#                 self.assertEqual(old_file,new_file)



if __name__ == '__main__':

    unittest.main()