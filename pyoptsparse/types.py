# Standard Python modules
from typing import Dict, List, Union

# External modules
from numpy import ndarray

NumpyType = Union[float, ndarray]
ArrayType = Union[float, List[float], ndarray]
Dict1DType = Dict[str, ndarray]
Dict2DType = Dict[str, Dict[str, ndarray]]
