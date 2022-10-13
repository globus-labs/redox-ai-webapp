"""Define the web application"""
from pathlib import Path
import logging

from fastapi import FastAPI
from fastapi.websockets import WebSocket, WebSocketDisconnect
from fastapi.staticfiles import StaticFiles
from fastapi.responses import FileResponse, HTMLResponse
from fastapi.middleware.cors import CORSMiddleware

from redoxweb.models import PropertyModel
from redoxweb.utils import make_svg_from_smiles
from redoxweb.config import models

app = FastAPI()
_my_dir = Path(__file__).parent
html_dir = _my_dir / 'html'
app.mount("/static", StaticFiles(directory=_my_dir / "static"), name="static")

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"]
)

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


@app.get('/api/models')
def list_models() -> list[PropertyModel]:
    """List out all models available with our service"""
    return models


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

            # Send a result back for each model
            for model in models:
                result = model.run(smiles)
                await socket.send_text(result.json())
    except WebSocketDisconnect:
        logger.info(f'Disconnected from client at {socket.client.host}')
