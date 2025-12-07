"""Legacy simulation module.

The functionality previously provided by ``SimulationRunner`` has been
superseded by the higher-level API in ``src.api``
(:func:`run_empirical_size_experiment`). This module is kept as a
placeholder to avoid import errors for any third-party code that might
still reference ``src.simulation``, but new code should use the API
module instead.
"""

__all__ = []
