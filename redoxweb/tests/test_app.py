from fastapi.testclient import TestClient

from redoxweb.app import app

test_app = TestClient(app)


def test_home():
    result = test_app.get("/")
    assert 'body' in result.text


def test_render():
    result = test_app.get("/api/render", params={'smiles': 'C'})
    assert 'svg' in result.text


def test_compute():
    with test_app.websocket_connect("/ws") as ws:
        ws.send_json({"smiles": "C"})
        data = ws.receive_json()
        assert data['smiles'] == "C"
