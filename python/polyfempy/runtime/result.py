from typing import Any, Mapping

import numpy as np
from numpy.typing import ArrayLike, NDArray


class Result:
    """The solution and metadata returned by a forward solve.

    ``sol`` is the raw solution returned by ``polyfem::State::solve``. It is
    not assumed to be aligned with mesh vertices or sampled visualization data.

    Attributes:
        sol: Contiguous NumPy array containing the raw solution coefficients.
            Its shape and component order depend on the finite-element space.
        meta: Dictionary of backend metadata. Its entries depend on the solve.
    """

    def __init__(self, sol: ArrayLike, *, meta: Mapping[str, Any] | None = None):
        self.sol: NDArray = np.ascontiguousarray(np.asarray(sol))
        self.meta: dict[str, Any] = {} if meta is None else dict(meta)
