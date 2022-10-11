"""Definition of the models to run"""

from rdkit import Chem
from rdkit.Chem import Crippen, Descriptors

from redoxweb.models import PropertyModel


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
