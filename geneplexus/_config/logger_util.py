"""Logger utilities."""
import logging
from typing import Optional

from .config import LOG_LEVEL_TYPE


def make_logger(
    name: Optional[str] = None,
    log_level: LOG_LEVEL_TYPE = "INFO",
) -> logging.Logger:
    """Make a basic logger used by GenePlexus.

    name (str, optional): Name of the logger.
    log_level (LOG_LEVEL): Level of details to log

    """
    formatter = logging.Formatter("%(name)s:%(funcName)s:%(levelname)s:%(message)s")
    stream_handler = logging.StreamHandler()
    stream_handler.setFormatter(formatter)
    stream_handler.setLevel(logging.DEBUG)

    logger_name = "geneplexus" if name is None else f"geneplexus.{name}"
    logger = logging.getLogger(logger_name)
    logger.addHandler(stream_handler)
    logger.setLevel(logging.getLevelName(log_level))

    return logger
