# Standard Python modules
from typing import Dict, Sequence, Union

# External modules
import numpy as np
import numpy.typing as npt

# Either ndarray or scalar
NumpyType = Union[float, npt.NDArray[np.float64]]
# ndarray, list of numbers, or scalar
ArrayType = Union[NumpyType, Sequence[float]]
# funcs
Dict1DType = Dict[str, npt.NDArray[np.float64]]
# funcsSens
Dict2DType = Dict[str, Dict[str, npt.NDArray[np.float64]]]
