import unittest
import numpy as np

# Force the local gqcpy to be imported
import sys
sys.path.insert(0, '../')

import gqcpy

class FukuiDysonAnalysisQCM(unittest.TestCase):

    def setUp(self):
        """ Iniates variables to be used by tests """
        N = gqcpy.Nucleus(7, 0, 0, 0)
        N_ = gqcpy.Nucleus(7, 2.1429494, 0, 0)
        N2 = gqcpy.Molecule([N,N_], +1) # The N2 molecule with an intramolecular distance of 2.1429494 bohr

        self.fukui_dyson_module = gqcpy.FukuiDysonAnalysis(N2, "STO-3G", False)

        self.fukui_vectors = np.array( [[-2.24456e-13, -3.28624e-14, -2.99164e-12,     0.994828,  9.90369e-12,     0.101575,   2.1121e-17,  4.71179e-17 , 3.27466e-14, -7.65573e-05],
                                        [0.0238145, -4.45007e-17, -3.48747e-16,  9.80093e-12,     -0.99971,  1.53623e-12,  3.09444e-12,  6.92411e-12,  -0.00360969,    -8.04077e-16],
                                        [6.44013e-14,  1.31846e-11,  3.40075e-10,    -0.101425,  5.32801e-13,     0.993397, -9.41798e-17,  2.07526e-16,  2.46607e-13,     0.0536178],
                                        [-0.511003, -1.18685e-16, -4.89684e-17,  2.24892e-14,   -0.0152761,  -1.3754e-13,  6.82643e-10,  1.53332e-09,     0.859443,     -5.9651e-13],
                                        [1.76243e-09,  6.25179e-12,  1.75868e-11,  5.47792e-17,  5.23169e-11, -1.01499e-15,     0.152336,     0.988329, -8.35439e-10,   1.54601e-14],
                                        [4.79937e-10, -1.76061e-11,  6.26495e-12, -3.54888e-17,  1.42587e-11,  2.11505e-16,     0.988329,    -0.152336, -2.27622e-10,  -1.63015e-15],
                                        [3.39825e-13, -1.23899e-10, -3.83883e-09,  -0.00552229,  3.84865e-14,    0.0533327,  7.66913e-16,  1.56215e-14, -4.81658e-13,     -0.998562],
                                        [5.9276e-16,     0.252821,     0.967513,  1.59139e-11, -4.00409e-16, -1.30116e-10, -4.42455e-12, -1.81337e-11,   4.6268e-16,   -3.75787e-09],
                                        [-3.44064e-17,     0.967513,    -0.252821, -3.44971e-12,   2.9481e-17,  2.72965e-11,  1.81566e-11, -4.41986e-12,  9.74385e-17,  8.53365e-10],
                                        [-0.859249,   2.1864e-16,  7.69967e-16, -1.13947e-14,   -0.0186227,  1.93388e-13,  4.58609e-10,  1.03042e-09,    -0.511218,    -3.61499e-14]] )

        self.fukui_eigenvalues = [-0.0225715, -0.00228606, -0.00228606, -8.02061e-07, 1.2942e-05, 0.000146312, 0.00812185, 0.00812185, 0.0746261, 0.936115]

        self.fukui_naturals = np.diagflat(self.fukui_eigenvalues)
        
        self.fukui_matrix = np.array(  [[7.21279e-07 -1.22858e-18  1.10019e-05 -4.49902e-16 -7.18738e-19  2.31245e-19  7.23604e-05  2.74003e-13 -6.24056e-14 -5.59942e-15],
                                        [-1.22858e-18  1.10586e-06  7.01173e-17  4.33624e-05 -6.63581e-13 -1.80578e-13  6.33228e-16  -1.1625e-19  6.73197e-20  0.000599824],
                                        [1.10019e-05  7.01173e-17   0.00283559  -1.3401e-14  7.77331e-16 -8.28171e-17   -0.0501125 -1.89396e-10  4.30038e-11  -9.9457e-15],
                                        [-4.49902e-16  4.33624e-05  -1.3401e-14    0.0492281 -2.01017e-11 -5.48084e-12  5.30627e-13  3.26706e-17  7.45734e-18   -0.0426986],
                                        [-7.18738e-19 -6.63581e-13  7.77331e-16 -2.01017e-11   0.00812185 -6.88171e-17 -1.43257e-14 -1.93545e-13 -1.66776e-14  7.48924e-11],
                                        [2.31245e-19 -1.80578e-13 -8.28171e-17 -5.48084e-12 -6.88171e-17   0.00812185  1.51052e-15 -1.67606e-14  1.93775e-13  2.03984e-11],
                                        [7.23604e-05  6.33228e-16   -0.0501125  5.30627e-13 -1.43257e-14  1.51052e-15     0.933425   3.5213e-09 -7.99644e-10  5.87618e-14],
                                        [2.74003e-13  -1.1625e-19 -1.89396e-10  3.26706e-17 -1.93545e-13 -1.67606e-14   3.5213e-09  -0.00228606 -4.46618e-17 -7.59975e-18],
                                        [-6.24056e-14  6.73197e-20  4.30038e-11  7.45734e-18 -1.66776e-14  1.93775e-13 -7.99644e-10 -4.46618e-17  -0.00228606 -4.30413e-18],
                                        [-5.59942e-15  0.000599824  -9.9457e-15   -0.0426986  7.48924e-11  2.03984e-11  5.87618e-14 -7.59975e-18 -4.30413e-18   0.00283837]] )

        self.dyson_coefficients = [ 9.71966e-05,  9.52394e-16,   -0.0516701,  6.91108e-13, -1.49184e-14,  1.57069e-15,    0.964144,  3.62447e-09, -8.23077e-10,  -4.4647e-14]                                

    def tearDown(self):
        pass

    def test_analysis(self):
        """ Compare the various analysis parameters with a reference value """
        self.assertAlmostEqual(self.fukui_dyson_module.get_fukui_naturals(),  self.fukui_naturals)
        self.assertAlmostEqual(self.fukui_dyson_module.get_fukui_matrix(),  self.fukui_matrix)
        self.assertAlmostEqual(self.fukui_dyson_module.get_dyson_coefficients(),  self.dyson_coefficients)
        self.assertAlmostEqual(self.fukui_dyson_module.get_fukui_vectors(),  self.fukui_vectors)

    
if __name__ == '__main__':
    unittest.main()
