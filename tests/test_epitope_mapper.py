import unittest
from chimerax.core.commands import run
from chimerax.epitope_mapper.engine import EpitopeScorer, SpatialClusterer

class TestEpitopeMapper(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        from chimerax.core import commands
        cls.session = commands.CommandOutput()

    def test_epitope_mapper(self):
        # Load a test model and map
        run(self.session, 'open C:\\Users\\Jerem\\Computational Physics\\ChimeraX-EpitopeMapper\\data\\testing\\Bxb1_Tail_Tip.cif')
        run(self.session, 'open C:\\Users\\Jerem\\Computational Physics\\ChimeraX-EpitopeMapper\\data\\testing\\emd_46661.map')

        # Select atoms
        run(self.session, 'select :1')

        # Run the epitope map command
        run(self.session, 'epitope map atoms :1 cutoff 5.0')

        # Add assertions to verify the output
        self.assertIn('Total residues evaluated', self.session.logger.messages)
        self.assertIn('Number of candidate residues exceeding the 1.5 threshold', self.session.logger.messages)
        self.assertIn('Region_1:', self.session.logger.messages)

if __name__ == '__main__':
    unittest.main()