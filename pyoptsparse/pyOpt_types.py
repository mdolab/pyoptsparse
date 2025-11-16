# Standard Python modules
from typing import Sequence

# External modules
import numpy as np
import numpy.typing as npt

# Either ndarray or scalar
NumpyType = float | npt.NDArray[np.float64]
# ndarray, list of numbers, or scalar
ArrayType = NumpyType | Sequence[float]
# funcs
Dict1DType = dict[str, npt.NDArray[np.float64]]
# funcsSens
Dict2DType = dict[str, dict[str, npt.NDArray[np.float64]]]
