"""Define the web application"""
from pathlib import Path
import logging

from fastapi import FastAPI
from fastapi.websockets import WebSocket, WebSocketDisconnect
from fastapi.staticfiles import StaticFiles
from fastapi.responses import FileResponse, HTMLResponse

from redoxweb.models import PredictionResult
from redoxweb.utils import make_svg_from_smiles

app = FastAPI()
_my_dir = Path(__file__).parent
html_dir = _my_dir / 'html'
app.mount("/static", StaticFiles(directory=_my_dir / "static"), name="static")

logger = logging.getLogger('app')


@app.get("/")
def read_root():
    """Serve the home page"""
    return FileResponse(html_dir / 'home.html')


@app.get('/api/render')
def render_molecule(smiles: str) -> HTMLResponse:
    """Render a molecule

    Args:
        smiles: SMILES string of the molecule to be rendered
    """
    return HTMLResponse(make_svg_from_smiles(smiles))


@app.websocket('/ws')
async def result_messenger(socket: WebSocket):
    """Open a socket connection that will manage sending results to-and-from the system

    Clients send the SMILES string of the molecule they want to evaluate and receive the results of
    computations as they finish.

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
