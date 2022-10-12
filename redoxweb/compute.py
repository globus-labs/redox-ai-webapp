"""Definition of the models to run"""
from functools import lru_cache
import os

import mdf_toolbox
from dlhub_sdk import DLHubClient
from rdkit import Chem
from rdkit.Chem import Crippen, Descriptors

from redoxweb.models import PropertyModel

_fx_scope = "https://auth.globus.org/scopes/facd7ccc-c5f4-42aa-916b-a0e270e2c2a9/all"


@lru_cache
def get_dlhub_client() -> DLHubClient:
    """Get a DLHub client

    Uses environment variables to create a confidential client, if available,
    and the normal token-based login procedure otherwise

    Returns:
        An already-authenticated client
    """

    if 'DLHUB_SECRET' in os.environ:
        # Pull the client information from the environment
        client_id = os.environ['DLHUB_CLIENT_ID']
        client_secret = os.environ['DLHUB_CLIENT_SECRET']

        # Get the services via a confidential log in
        services = ["search", "dlhub", _fx_scope, "openid"]
        auth_res = mdf_toolbox.confidential_login(client_id=client_id,
                                                  client_secret=client_secret,
                                                  services=services,
                                                  make_clients=False)
        return DLHubClient(
            dlh_authorizer=auth_res["dlhub"], fx_authorizer=auth_res[_fx_scope],
            openid_authorizer=auth_res['openid'], search_authorizer=auth_res['search'],
            force_login=False, http_timeout=10
        )
    else:
        return DLHubClient()


class LogPModel(PropertyModel):
    """Compute LogP using Crippen's method"""

    id = "logp"
    units = ""
    name = "LogP"

    def _run(self, smiles: str) -> float:
        mol = Chem.MolFromSmiles(smiles)
        return Crippen.MolLogP(mol)


class MolecularWeightModel(PropertyModel):
    """Compute the molecular weight of a molecule"""

    id = "molwt"
    units = "g/mol"
    name = "Molecular Weight"

    def _run(self, smiles: str) -> float:
        mol = Chem.MolFromSmiles(smiles)
        return Descriptors.MolWt(mol)


class SolvationEnergyModel(PropertyModel):
    """Run the solvation energy model from our 2021 paper"""
    id = "solv_ml"
    units = "kcal/mol"
    name = "G<sub>solv,ACN</sub>"

    def _run(self, smiles: str) -> float:
        client = get_dlhub_client()
        output = client.run("loganw_globusid/mpnn-solvation-energy", inputs=[smiles])
        return output['solvation-energies'][0][1]  # ACN is the second solvent
