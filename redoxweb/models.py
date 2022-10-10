from typing import Optional, Tuple

from pydantic import BaseModel, Field


class PredictionResult(BaseModel):
    """Result of inference on a DLHub model"""

    smiles: str = Field(..., help='SMILES string of molecule passed to the client')
    key: str = Field(..., help='InChI key of the molecule')
    model_name: str = Field(..., help='Name of the model which was invoked')
    value: float = Field(..., help='Output value')
    confidence_interval: Optional[Tuple[float, float]] = Field(None, help='Confidence intervals, if known')
