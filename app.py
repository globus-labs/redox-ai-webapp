from typing import Tuple, Optional
from pathlib import Path
import logging

from pydantic import BaseModel, Field
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import MolDraw2DSVG
from fastapi import FastAPI
from fastapi.websockets import WebSocket, WebSocketDisconnect
from fastapi.staticfiles import StaticFiles
from fastapi.responses import FileResponse, HTMLResponse

app = FastAPI()
html_dir = Path(__file__).parent / 'html'
app.mount("/static", StaticFiles(directory="static"), name="static")

logger = logging.getLogger('app')


class PredictionResult(BaseModel):
    """Result of inference on a DLHub model"""

    smiles: str = Field(..., help='SMILES string of molecule passed to the client')
    key: str = Field(..., help='InChI key of the molecule')
    model_name: str = Field(..., help='Name of the model which was invoked')
    value: float = Field(..., help='Output value')
    confidence_interval: Optional[Tuple[float, float]] = Field(None, help='Confidence intervals, if known')


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
    """Serve the home page"""
    return FileResponse(html_dir / 'home.html')


@app.get('/api/render')
def render_molecule(smiles: str) -> HTMLResponse:
    """Render a molecule"""
    return HTMLResponse(_print_molecule(smiles))


@app.websocket('/ws')
async def result_messenger(socket: WebSocket):
    """Open a socket connection that will manage sending results to-and-from the system

    Args:
        socket: The websocket created for this particular session
    """

    # Accept the connection
    await socket.accept()
    logger.info(f'Connected to client at {socket.client.host}')

    try:
        while True:
            # Wait for messages from the webapp to arrive
            msg = await socket.receive_json()
            smiles = msg['smiles']

            # Send a result back
            result = PredictionResult(
                smiles=smiles,
                key='aaa' * 4,
                model_name='solv_ml',
                value=1.
            )
            await socket.send_text(result.json())
    except WebSocketDisconnect:
        logger.info(f'Disconnected from client at {socket.client.host}')
