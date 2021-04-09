# Standard Python modules
from typing import Dict, List, Union

# External modules
from numpy import ndarray

# Either ndarray or scalar
NumpyType = Union[float, ndarray]
# ndarray, list of numbers, or scalar
ArrayType = Union[float, List[float], ndarray]
# funcs
Dict1DType = Dict[str, ndarray]
# funcsSens
Dict2DType = Dict[str, Dict[str, ndarray]]
