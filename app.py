from pathlib import Path
import logging

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import MolDraw2DSVG
from fastapi import FastAPI
from fastapi.staticfiles import StaticFiles
from fastapi.responses import FileResponse, HTMLResponse

app = FastAPI()
html_dir = Path(__file__).parent / 'html'
app.mount("/static", StaticFiles(directory="static"), name="static")

logger = logging.getLogger('app')


def _print_molecule(smiles) -> str:
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


@app.get("/")
def read_root():
    return FileResponse(html_dir / 'home.html')


@app.get('/api/render')
def render_molecule(smiles: str):
    return HTMLResponse(_print_molecule(smiles))
