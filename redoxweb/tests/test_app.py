from fastapi.testclient import TestClient
from pytest import mark

from redoxweb.app import app
from redoxweb.config import models

test_app = TestClient(app)


def test_home():
    result = test_app.get("/")
    assert 'body' in result.text


def test_render():
    result = test_app.get("/api/render", params={'smiles': 'C'})
    assert 'svg' in result.text


def test_models():
    result = test_app.get("/api/models")
    assert len(result.json()) == len(models)


@mark.timeout(10)
def test_compute():
    with test_app.websocket_connect("/ws") as ws:
        ws.send_json({"smiles": "C"})

        # Make sure it sends a result back for each model
        seen = set()
        known = set(m.id for m in models)
        for _ in models:
            data = ws.receive_json()
            assert data['smiles'] == "C"
            seen.add(data['model_name'])
        assert seen == known
