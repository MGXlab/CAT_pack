import errno
import os
from contextlib import contextmanager, ExitStack
from pathlib import Path
from typing import Iterable, Iterator

from .errors import InputError


@contextmanager
def _lock_windows(path: Path) -> Iterator[None]:
    # O_TEMPORARY: Windows will delete the lockfile itself once the last handle has been closed
    flags = os.O_RDWR | os.O_CREAT | os.O_TEMPORARY | os.O_BINARY

    # DEVNOTE: python versions 3.15, 3.12, 3.14 do not have module msvcrt on Linux, so only import after _lock_windows has been selected
    import msvcrt

    file = os.fdopen(os.open(path, flags, 0o600), "r+b")
    try:
        msvcrt.locking(file.fileno(), msvcrt.LK_NBLCK, 1)
        try:
            yield
        finally:
            file.seek(0)
            msvcrt.locking(file.fileno(), msvcrt.LK_UNLCK, 1)
    finally:
        file.close()


@contextmanager
def _lock_posix(path: Path) -> Iterator[None]:
    import fcntl

    flags = os.O_RDWR | os.O_CREAT

    # DEVNOTE: Due to differences in file locking between windows and linux
    #   we have to retry if the lock file got deleted (or replaced) when the current
    #   process is trying to place a lock. (race condition)
    while True:
        file = os.fdopen(os.open(path, flags, 0o600), "r+b")
        try:
            fcntl.flock(file.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
            try:
                try:
                    same_file = os.path.samestat(os.fstat(file.fileno()), path.stat())
                except FileNotFoundError:
                    same_file = False

                if not same_file:
                    # If a previous process deleted it, the file will not show up as a same file
                    # Thus we must try again
                    continue

                try:
                    yield
                finally:
                    # Delete while the output is still locked so no other process grab edit/remove the "old" file
                    path.unlink(missing_ok=True)
            finally:
                fcntl.flock(file.fileno(), fcntl.LOCK_UN)
        finally:
            file.close()
        return


@contextmanager
def _lock_file(path: Path) -> Iterator[None]:
    if os.name == "nt":
        with _lock_windows(path):
            yield
    else:
        with _lock_posix(path):
            yield

@contextmanager
def lock_outputs(resources: Iterable[Path]) -> Iterator[None]:
    """Lock dedicated lock files and remove them when their ownership ends."""
    paths = list(dict.fromkeys(p.resolve() for p in resources))

    with ExitStack() as stack:
        for path in paths:
            try:
                stack.enter_context(_lock_file(path))
            except OSError as error:
                if error.errno in (errno.EACCES, errno.EAGAIN, errno.EDEADLK):
                    raise InputError(
                        "Another run is using the output location or prefix.", path=path,
                        hint="Wait for that run to finish or choose a different output location.",
                    )
                else:
                    raise InputError("Could not lock the output location.", path=path, hint="CAT pack couldn't reserve this output location. Check that you can write to the folder, the folder isn't used by another process, or choose a different output folder/location")
        yield
