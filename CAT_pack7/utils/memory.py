import logging

import psutil

from .errors import InputError
from ..config.settings import MemoryBudget, MemorySettings

log = logging.getLogger("CAT_pack")
MIB = 1024 ** 2
GIB = 1024 ** 3
PROCESS_RESERVE = 128 * MIB




def get_memory_budget(settings: MemorySettings, workers: int) -> MemoryBudget:
    available = psutil.virtual_memory().available
    requested = settings.available_memory_bytes
    if requested is None:
        requested = min(GIB, available) if settings.low_memory else available
    total = min(requested, available)
    if requested > available:
        log.warning("Requested memory budget exceeds detected available RAM; using %.1f MiB.", total / MIB)
    parent_reserve = PROCESS_RESERVE
    workers = min(workers, (total - parent_reserve - MIB) // PROCESS_RESERVE)
    if workers < 1:
        raise InputError(
            "Too little available memory for the estimated process overhead.",
            hint="Free memory or increase --available-memory; allow at least 128 MiB each for parent and worker plus sequence storage.",
        )
    batch_bytes = total - parent_reserve - workers * PROCESS_RESERVE
    return MemoryBudget(total, workers, batch_bytes)
