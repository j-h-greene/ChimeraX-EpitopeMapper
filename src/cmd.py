from chimerax.core.commands import CmdDesc, AtomsArg, FloatArg
from chimerax.map import volume_from_grid_data
from chimerax.epitope_mapper.engine import EpitopeScorer, SpatialClusterer
from chimerax.core.models import Atom
from chimerax.map_data import sample_volume_data
from chimerax.core.errors import UserError
import numpy as np

def epitope_map(session, atoms, cutoff=5.0):
    try:
        # Check if the selection contains protein atoms
        if not any(isinstance(atom, Atom) for atom in atoms):
            raise UserError("The selection must contain protein atoms.")

        # Check if there is a density map loaded
        if not session.maps:
            raise UserError("No density map data is currently loaded. Please load a density map to sample from.")

        # Extract unique residues from the passed atoms
        unique_residues = {atom.residue for atom in atoms}

        residues_data = []
        for residue in unique_residues:
            amino_acid_type = residue.name
            chain = residue.chain_id
            residue_id = residue.number
            coordinates = residue.atoms.principal_axes().center

            # Calculate the average map value for the residue
            map_values = []
            for atom in residue.atoms:
                map_value = sample_volume_data(session.maps[0], atom.coord)
                map_values.append(map_value)
            raw_map_value = np.mean(map_values)

            residues_data.append({
                'residue_id': residue_id,
                'chain': chain,
                'amino_acid_type': amino_acid_type,
                'raw_map_value': raw_map_value,
                'coordinates': coordinates
            })

        # Feed the extracted dataset into EpitopeScorer and SpatialClusterer
        scorer = EpitopeScorer(residues_data)
        high_scoring_residues = scorer.get_sorted_residues_above_threshold(threshold=1.5)
        clusterer = SpatialClusterer(high_scoring_residues)
        regions = clusterer.find_conformational_regions(distance_cutoff=cutoff)

        # Log the results
        session.logger.info(f"Total residues evaluated: {len(residues_data)}")
        session.logger.info(f"Number of candidate residues exceeding the 1.5 threshold: {len(high_scoring_residues)}")

        # Sort regions by Total Epitope Footprint descending
        sorted_regions = sorted(regions.items(), key=lambda item: item[1]['total_epitope_footprint'], reverse=True)

        for region_id, metrics in sorted_regions:
            session.logger.info(f"<b>{region_id}:</b>")
            session.logger.info(f"  Average Region Strength: {metrics['average_region_strength']:.2f}")
            session.logger.info(f"  Total Epitope Footprint: {metrics['total_epitope_footprint']:.2f}")
            session.logger.info("  Residues:")
            for residue in metrics['residues']:
                session.logger.info(f"    Chain {residue['chain']}: {residue['amino_acid_type']} {residue['residue_id']} - Epitope Score: {residue['epitope_score']:.2f}")

    except UserError as e:
        session.logger.error(str(e))

epitope_map_desc = CmdDesc(
    required=[
        ('atoms', AtomsArg),
    ],
    keyword=[
        ('cutoff', FloatArg(default=5.0)),
    ],
    synopsis='Map epitopes on protein structures.'
)