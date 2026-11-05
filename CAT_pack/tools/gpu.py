import logging
import sys
from pathlib import Path

from ..utils.errors import ExternalToolError, InputError
from ..utils.memory import GIB

log = logging.getLogger("CAT_pack")


def check_mmseqs_gpu(database: Path) -> None:
    # DEVNOTE: The reason for linux only support is because checking on windows
    #  For CUDA has proven to be quite difficult
    if sys.platform != "linux":
        raise InputError("MMseqs2 GPU search is only supported on Linux")
    try:
        import pycuda.driver as cuda
    except ImportError as error:
        raise ExternalToolError("CUDA", f"PyCUDA is unavailable: {error}",
                                hint="Please Install CAT_pack's GPU extra: pip install '.[gpu]'.")

    try:
        cuda.init()
        if cuda.Device.count() == 0:
            raise InputError("No CUDA GPUs are available.")
        for ordinal in range(cuda.Device.count()):
            device = cuda.Device(ordinal)
            if device.compute_capability() < (7, 5):
                raise InputError(f"GPU {ordinal} ({device.name()}) is too old for MMseqs2.",
                                 hint="See MMseqs2 documentation for supported GPU's")
            context = device.make_context()
            try:
                free, total = cuda.mem_get_info()
            finally:
                try:
                    context.pop()
                finally:
                    context.detach()
            log.info("GPU %d (%s): %.2f GiB free / %.2f GiB total.",
                     ordinal, device.name(), free / GIB, total / GIB)
    except cuda.Error as error:
        raise ExternalToolError("CUDA", f"GPU availability check failed: {error}")


