"""Misc utility functions"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import MolDraw2DSVG


def make_svg_from_smiles(smiles) -> str:
    """Print a molecule as an SVG

    Args:
        smiles (str): SMILES string of molecule to present
    Returns:
        (str): SVG rendering of molecule
    """

    # Compute 2D coordinates
    mol = Chem.MolFromSmiles(smiles)
    AllChem.Compute2DCoords(mol)

    # Print out an SVG
    rsvg = MolDraw2DSVG(256, 256)
    rsvg.DrawMolecule(mol)
    rsvg.FinishDrawing()
    return rsvg.GetDrawingText().strip()
