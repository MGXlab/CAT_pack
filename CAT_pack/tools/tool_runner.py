"""Run external tools

Run external tools with the error and log handling

Here for possible future buildout with an extra layer of log handling
"""
import logging
import os
import subprocess
from contextlib import contextmanager
from io import BufferedWriter
from threading import Thread
from typing import Generator

from ..utils.errors import ExternalToolError

log = logging.getLogger("CAT_pack")


@contextmanager
def tool_output(tool: str) -> Generator[BufferedWriter, None, None]:
    read_fd, write_fd = os.pipe()
    reader = os.fdopen(read_fd, "rb")
    writer = os.fdopen(write_fd, "wb")

    def drain():
        with reader:
            for line in reader:
                log.info("[%s] %s", tool, line.decode(errors="replace").rstrip())

    thread = Thread(target=drain, name=f"{tool}-output", daemon=True)
    thread.start()
    try:
        yield writer
    finally:
        writer.close()
        thread.join()


def run_tool(command: list[str], tool: str, operation: str = "command") -> None:
    log.info("Running command: %s", " ".join(command))
    try:
        with tool_output(tool) as output:
            subprocess.run(command, check=True, stdout=output,
                           stderr=subprocess.STDOUT)
    except (OSError, subprocess.CalledProcessError) as error:
        raise ExternalToolError(tool, f"{operation} failed: {error}") from error
