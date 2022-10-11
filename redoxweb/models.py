"""Data models describing inputs and outputs of the web application"""

from typing import Optional, Tuple

from pydantic import BaseModel, Field


class PredictionResult(BaseModel):
    """Result of inference on a DLHub model"""

    smiles: str = Field(..., help='SMILES string of molecule passed to the client')
    model_name: str = Field(..., help='Name of the model which was invoked')
    value: float = Field(..., help='Output value')
    confidence_interval: Optional[Tuple[float, float]] = Field(None, help='Confidence intervals, if known')


class PropertyModel(BaseModel):
    """Base class for a machine learning model that predicts a molecular property from SMILES"""

    id: str = Field(..., help='Short name of the molecules')
    name: str = Field(..., help='Name of the property being evaluated')
    units: str = Field(..., help='Units of the molecule')

    def _run(self, smiles: str) -> float:
        """Predict the property of a molecule

        Args:
            smiles: SMILES of the molecule being evaluated
        Returns:
            Predicted property
        """
        raise NotImplementedError()

    def run(self, smiles: str) -> PredictionResult:
        """Run the model and return a formatted result object"""

        return PredictionResult(
            smiles=smiles,
            model_name=self.id,
            value=self._run(smiles)
        )
