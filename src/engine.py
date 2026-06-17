import numpy as np
from scipy.ndimage import gaussian_filter
from scipy.spatial import KDTree
from collections import defaultdict

class EpitopeScorer:
    def __init__(self, residues):
        self.residues = residues
        self.gaussian_map_values = []
        self.protrusion_scores = []
        self.propensity_scores = []
        self.epitope_scores = []

    def _calculate_gaussian_map_values(self):
        raw_values = [residue['raw_map_value'] for residue in self.residues]
        self.gaussian_map_values = gaussian_filter(raw_values, sigma=1)

    def _normalize_gaussian_map_values(self):
        max_value = max(self.gaussian_map_values)
        min_value = min(self.gaussian_map_values)
        if max_value == min_value:
            self.protrusion_scores = [0.5] * len(self.gaussian_map_values)
        else:
            self.protrusion_scores = [(max_value - value) / (max_value - min_value) for value in self.gaussian_map_values]

    def _calculate_propensity_scores(self):
        propensity_scale = {
            'TYR': 1.0, 'TRP': 1.0, 'PHE': 1.0,
            'ARG': 0.9, 'LYS': 0.9, 'HIS': 0.9,
            'ASP': 0.8, 'GLU': 0.8, 'ASN': 0.8, 'GLN': 0.8,
            'SER': 0.6, 'THR': 0.6, 'GLY': 0.6, 'PRO': 0.6,
            'CYS': 0.3, 'MET': 0.3,
            'ALA': 0.1, 'VAL': 0.1, 'LEU': 0.1, 'ILE': 0.1
        }
        self.propensity_scores = [propensity_scale.get(residue['amino_acid_type'], 0.0) for residue in self.residues]

    def _calculate_epitope_scores(self):
        self.epitope_scores = [p + g for p, g in zip(self.propensity_scores, self.protrusion_scores)]

    def calculate_scores(self):
        self._calculate_gaussian_map_values()
        self._normalize_gaussian_map_values()
        self._calculate_propensity_scores()
        self._calculate_epitope_scores()

    def get_sorted_residues_above_threshold(self, threshold=1.5):
        self.calculate_scores()
        scored_residues = [
            {**residue, 'epitope_score': score}
            for residue, score in zip(self.residues, self.epitope_scores)
        ]
        return sorted(
            [residue for residue in scored_residues if residue['epitope_score'] > threshold],
            key=lambda x: x['epitope_score'],
            reverse=True
        )

class SpatialClusterer:
    def __init__(self, residues):
        self.residues = residues
        self.coordinates = [residue['coordinates'] for residue in residues]
        self.kdtree = KDTree(self.coordinates)

    def find_conformational_regions(self, distance_cutoff=5.0):
        regions = defaultdict(list)
        visited = set()

        def find_neighbors(index):
            neighbors = self.kdtree.query_ball_point(self.coordinates[index], distance_cutoff)
            return neighbors

        region_id = 1
        for i in range(len(self.residues)):
            if i in visited:
                continue
            stack = [i]
            while stack:
                current = stack.pop()
                if current in visited:
                    continue
                visited.add(current)
                regions[region_id].append(self.residues[current])
                neighbors = find_neighbors(current)
                for neighbor in neighbors:
                    if neighbor not in visited:
                        stack.append(neighbor)
            region_id += 1

        return self._calculate_region_metrics(regions)

    def _calculate_region_metrics(self, regions):
        region_metrics = {}
        for region_id, residues in regions.items():
            epitope_scores = [residue['epitope_score'] for residue in residues]
            average_strength = np.mean(epitope_scores)
            total_footprint = np.sum(epitope_scores)
            region_metrics[f'Region_{region_id}'] = {
                'residues': residues,
                'average_region_strength': average_strength,
                'total_epitope_footprint': total_footprint
            }
        return region_metrics

if __name__ == "__main__":
    import unittest

    class TestEpitopeScorer(unittest.TestCase):
        def setUp(self):
            self.residues = [
                {'residue_id': 1, 'chain': 'A', 'amino_acid_type': 'TYR', 'raw_map_value': 0.5, 'coordinates': (0.0, 0.0, 0.0)},
                {'residue_id': 2, 'chain': 'A', 'amino_acid_type': 'ARG', 'raw_map_value': 0.3, 'coordinates': (1.0, 1.0, 1.0)},
                {'residue_id': 3, 'chain': 'A', 'amino_acid_type': 'ASP', 'raw_map_value': 0.7, 'coordinates': (2.0, 2.0, 2.0)},
                {'residue_id': 4, 'chain': 'A', 'amino_acid_type': 'CYS', 'raw_map_value': 0.9, 'coordinates': (3.0, 3.0, 3.0)},
                {'residue_id': 5, 'chain': 'A', 'amino_acid_type': 'ALA', 'raw_map_value': 0.1, 'coordinates': (4.0, 4.0, 4.0)},
            ]

        def test_calculate_scores(self):
            scorer = EpitopeScorer(self.residues)
            scorer.calculate_scores()
            self.assertEqual(len(scorer.epitope_scores), len(self.residues))
            self.assertAlmostEqual(sum(scorer.epitope_scores), 5.9, places=1)

        def test_get_sorted_residues_above_threshold(self):
            scorer = EpitopeScorer(self.residues)
            sorted_residues = scorer.get_sorted_residues_above_threshold(threshold=1.5)
            self.assertEqual(len(sorted_residues), 2)
            self.assertEqual(sorted_residues[0]['residue_id'], 1)
            self.assertEqual(sorted_residues[1]['residue_id'], 3)

        def test_spatial_clustering(self):
            scorer = EpitopeScorer(self.residues)
            high_scoring_residues = scorer.get_sorted_residues_above_threshold(threshold=1.5)
            clusterer = SpatialClusterer(high_scoring_residues)
            regions = clusterer.find_conformational_regions(distance_cutoff=2.0)
            self.assertEqual(len(regions), 2)
            self.assertIn('Region_1', regions)
            self.assertIn('Region_2', regions)
            self.assertEqual(len(regions['Region_1']['residues']), 1)
            self.assertEqual(len(regions['Region_2']['residues']), 1)

    unittest.main()