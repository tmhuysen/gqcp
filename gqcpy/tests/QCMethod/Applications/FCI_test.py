import unittest
import numpy as np

# Force the local gqcpy to be imported
import sys
sys.path.insert(0, '../')

import gqcpy


class FCITest(unittest.TestCase):

    def setUp(self):
        """Initiate variables to be used by the tests"""

        # Provide the reference data for H2O
        self.fci_module = gqcpy.FCI("data/h2o_Psi4_GAMESS.xyz", "STO-3G", 5, 5, False)
        self.ref_energy = -75.0129803939602
        self.fci_module.solve()


    def tearDown(self):
        pass


    def test_energies(self):
        """Compare the energy with a reference value"""
        self.assertAlmostEqual(self.fci_module.get_energy(), self.ref_energy)


if __name__ == '__main__':
    unittest.main()
