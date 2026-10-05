"""Run external tools

Run external tools with the error and log handling

Here for possible future buildout with an extra layer of log handling
"""
import logging
import subprocess

from ..utils.errors import ExternalToolError

log = logging.getLogger("CAT_pack")


def run_tool(command: list[str], tool: str, operation: str = "command") -> None:
    log.info("Running command: %s", " ".join(command))
    try:
        subprocess.run(command, check=True)
    except (OSError, subprocess.CalledProcessError) as error:
        raise ExternalToolError(tool, f"{operation} failed: {error}") from error
