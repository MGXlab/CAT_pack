import errno
import os
from contextlib import contextmanager, ExitStack
from pathlib import Path
from typing import Iterable, Iterator

from .errors import InputError


@contextmanager
def lock_outputs(resources: Iterable[Path]) -> Iterator[None]:
    with ExitStack() as stack:
        for path in dict.fromkeys(path.resolve() for path in resources):
            lock_path = path.with_name(f".{path.name}.from.cat.lock")
            lock_file = stack.enter_context(lock_path.open("a+b"))
            try:
                lock_file.seek(0)
                if os.name == "nt": # <- For future reference: THis is for windows
                    import msvcrt
                    msvcrt.locking(lock_file.fileno(), msvcrt.LK_NBLCK, 1)
                else:
                    import fcntl
                    fcntl.flock(lock_file.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
            except OSError as error:
                if error.errno in {errno.EACCES, errno.EAGAIN, errno.EDEADLK}:
                    raise InputError(
                        "Another run is using the output location or prefix.", path=path,
                        hint="Wait for that run to finish or choose a different output location.",
                    )
                else:
                    raise InputError("Could not lock the output location.", path=path)
        yield
