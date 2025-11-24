from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Optional


@dataclass
class RuntimeOptions:
    device: str = "cpu"  # 'cpu' | 'cuda'
    dtype: str = "float32"
    chunk_size: Optional[int] = None
    threads: Optional[int] = None
    seed: Optional[int] = None


def options_from_basicio(basicio: Any) -> RuntimeOptions:
    return RuntimeOptions(
        device=getattr(basicio, "device", "cpu"),
        dtype=getattr(basicio, "dtype", "float32"),
        chunk_size=getattr(basicio, "chunk_size", None),
        threads=getattr(basicio, "threads", None),
        seed=getattr(basicio, "seed", None),
    )


