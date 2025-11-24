from __future__ import annotations

from typing import Any, Dict


def write_outputs_placeholder(GalaxyIO: Any) -> Dict[str, Any]:
    """Placeholder to keep MakeCubes flow compatible while we implement.

    Returns a dictionary of file paths/metadata that higher layers might expect.
    """
    return {"output_folder": getattr(GalaxyIO, "OutputFolder", None)}


