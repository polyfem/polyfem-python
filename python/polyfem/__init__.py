"""PolyFEM model builder and generated configuration factories."""

from typing import TYPE_CHECKING

from polyfempy.generated_api import generated_api as _api
from polyfempy.generated_api.generated_api import Root as Root
from polyfempy.generated_api.generated_api import model as model

if TYPE_CHECKING:
    from polyfempy.generated_api.generated_api import *  # noqa: F403
else:
    __all__ = ["Root", *(name for name, value in vars(_api).items()
                         if not name.startswith("_") and getattr(value, "__module__", None) == _api.__name__)]
    globals().update({name: getattr(_api, name) for name in __all__})
