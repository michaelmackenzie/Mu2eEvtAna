"""Shared Mu2e constants from the Mu2eEvtAna plotter's ``physics/Mu2ePhysics.C`` (the single source of these numbers).

    from mu2ecombine import mu2e_physics as phys
    npot = phys.run1a.npot
    f_cap = phys.muon_capture_fraction_al

The header is loaded into ROOT on first use; each attribute access returns a float, or a nested namespace.
"""

import os


class PhysicsError(AttributeError):
    pass


def physics_header():
    here = os.path.dirname(os.path.abspath(__file__))
    return os.path.normpath(os.path.join(here, "..", "..", "plotter", "physics", "Mu2ePhysics.C"))


_loaded = False


def _root_namespace():
    global _loaded
    import ROOT
    if not _loaded:
        header = physics_header()
        if not os.path.exists(header):
            raise PhysicsError(f"Mu2e physics constants not found at {header}")
        if not ROOT.gInterpreter.Declare(f'#include "{header}"'):
            raise PhysicsError(f"Unable to load {header}")
        _loaded = True
    return ROOT.mu2e_physics


class _Namespace:
    def __init__(self, path=()):
        self._path = path

    def __getattr__(self, name):
        if name.startswith("_"):
            raise AttributeError(name)
        obj = _root_namespace()
        path = self._path + (name,)
        try:
            for part in path:
                obj = getattr(obj, part)
        except AttributeError:
            raise PhysicsError(f"Unknown constant mu2e_physics.{'.'.join(path)}") from None
        if isinstance(obj, (int, float)):
            return float(obj)
        try:
            return float(obj)
        except (TypeError, ValueError):
            pass
        if callable(obj):
            return obj
        return _Namespace(path)

    def __repr__(self):
        return f"mu2e_physics{'.' if self._path else ''}{'.'.join(self._path)}"


mu2e_physics = _Namespace()
