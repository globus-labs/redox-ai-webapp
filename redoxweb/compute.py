"""Definition of the models to run"""

from rdkit import Chem
from rdkit.Chem import Crippen, Descriptors

from redoxweb.models import PropertyModel


class LogPModel(PropertyModel):
    """Compute LogP using Crippen's method"""

    name = "logp"
    units = ""

    def _run(self, smiles: str) -> float:
        mol = Chem.MolFromSmiles(smiles)
        return Crippen.MolLogP(mol)


class MolecularWeightModel(PropertyModel):
    """Compute the molecular weight of a molecule"""

    name = "molwt"
    units = "g/mol"

    def _run(self, smiles: str) -> float:
        mol = Chem.MolFromSmiles(smiles)
        return Descriptors.MolWt(mol)
