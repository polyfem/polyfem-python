import importlib
import numpy as np


def _maybe_import(name: str):
    """Import a module only when needed.

    Args:
        name: Module name, e.g. "torch" or "jax.numpy".

    Returns:
        The imported module object if available; otherwise None.
    """
    try:
        return importlib.import_module(name)
    except Exception:
        return None


def _is_torch_tensor(x) -> bool:
    """Check whether `x` is a torch.Tensor without importing torch.

    We detect by inspecting the object's type module name.

    Args:
        x: Any Python object.

    Returns:
        True if `x` looks like a torch.Tensor; False otherwise.
    """
    m = type(x).__module__
    return m.startswith("torch") and type(x).__name__ == "Tensor"


def _is_jax_array(x) -> bool:
    """Check whether `x` is a JAX array (no import required).

    Args:
        x: Any Python object.

    Returns:
        True if `x` originates from a 'jax.*' module; False otherwise.
    """
    m = type(x).__module__
    return m.startswith("jax")


def detect_backend(x):
    """Detect the backend of `x` without forcing imports.

    Args:
        x: Array/tensor-like object from NumPy, Torch, or JAX.

    Returns:
        Literal string: 'numpy', 'torch', or 'jax'.
    """
    if _is_torch_tensor(x):
        return "torch"
    if _is_jax_array(x):
        return "jax"
    return "numpy"


def _torch_to_numpy_zero_copy(t):
    """Convert a Torch tensor to NumPy with zero copy when safe.

    Conditions for zero copy:
    - CPU device
    - contiguous memory layout

    Falls back to moving to CPU and making it contiguous if needed.

    Args:
        t: torch.Tensor (any device/contiguity).

    Returns:
        np.ndarray sharing memory with `t` if (CPU + contiguous), otherwise a CPU copy.

    Raises:
        RuntimeError: If a torch tensor is detected but torch is not installed.
    """
    torch = _maybe_import("torch")
    if torch is None:
        raise RuntimeError("Torch tensor detected but 'torch' is not installed.")
    t = t.detach()
    if t.device.type != "cpu":
        t = t.cpu()
    if not t.is_contiguous():
        t = t.contiguous()
    return t.numpy()  # shares memory with `t` when contiguous CPU


def as_numpy(x, dtype=None):
    """Normalize any backend array/tensor to a C-contiguous NumPy array on CPU.

    Torch CPU contiguous tensors go through zero copy; Torch GPU tensors are moved
    to CPU; JAX arrays are converted via np.asarray (may copy).

    Args:
        x: Input array/tensor (NumPy / Torch / JAX).
        dtype: Optional target dtype; if provided, cast without copying when possible.

    Returns:
        (arr_np, backend):
            arr_np (np.ndarray): CPU, C-contiguous NumPy array.
            backend (str): 'numpy' | 'torch' | 'jax' indicating the original backend.
    """
    backend = detect_backend(x)
    if backend == "torch":
        arr = _torch_to_numpy_zero_copy(x)
    elif backend == "jax":
        arr = np.asarray(x)
    else:
        arr = np.asarray(x)

    if dtype is not None and arr.dtype != dtype:
        arr = arr.astype(dtype, copy=False)
    if not arr.flags.c_contiguous:
        arr = np.ascontiguousarray(arr)
    return arr, backend


def from_numpy(arr, backend):
    """Convert a NumPy array back to the requested backend.

    Torch path uses zero copy (shared memory) when possible.

    Args:
        arr: NumPy array (ideally C-contiguous for torch zero copy).
        backend: Target backend: 'numpy' | 'torch' | 'jax'.

    Returns:
        Array/tensor in the target backend.

    Raises:
        RuntimeError: If torch/jax is requested but not installed.
    """
    if backend == "torch":
        torch = _maybe_import("torch")
        if torch is None:
            raise RuntimeError("Requested torch backend, but 'torch' is not installed.")
        if not arr.flags.c_contiguous:
            arr = np.ascontiguousarray(arr)
        return torch.from_numpy(arr) 

    if backend == "jax":
        jnp = _maybe_import("jax.numpy")
        if jnp is None:
            raise RuntimeError("Requested jax backend, but 'jax' is not installed.")
        return jnp.asarray(arr)

    return arr  


def to_backend(arr, backend):
    """None-safe wrapper around from_numpy().

    Args:
        arr: NumPy array or None (when a field/result is optional).
        backend: Target backend for conversion.

    Returns:
        None if `arr` is None; otherwise converted array/tensor.
    """
    if arr is None:
        return None
    return from_numpy(arr, backend)
  