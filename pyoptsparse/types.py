# Standard Python modules
from typing import Dict, List, Union

# External modules
import numpy as np
import numpy.typing as npt

# Either ndarray or scalar
NumpyType = Union[float, npt.NDArray[np.float_]]
# ndarray, list of numbers, or scalar
ArrayType = Union[NumpyType, List[float]]
# funcs
Dict1DType = Dict[str, npt.NDArray[np.float_]]
# funcsSens
Dict2DType = Dict[str, Dict[str, npt.NDArray[np.float_]]]
