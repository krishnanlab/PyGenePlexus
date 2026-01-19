"""Logger utilities."""

import logging
from contextlib import contextmanager
from typing import Optional

from .config import LOG_LEVEL_TYPE


def make_logger(
    name: Optional[str] = None,
    log_level: LOG_LEVEL_TYPE = "INFO",
) -> logging.Logger:
    """Make a basic logger used by GenePlexus.

    name: Name of the logger.
    log_level: Level of details to log

    """
    formatter = logging.Formatter("%(name)s:%(funcName)s:%(levelname)s:%(message)s")
    stream_handler = logging.StreamHandler()
    stream_handler.setFormatter(formatter)
    stream_handler.setLevel(logging.getLevelName(log_level))

    logger_name = "geneplexus" if name is None else f"geneplexus.{name}"
    logger = logging.getLogger(logger_name)
    logger.addHandler(stream_handler)
    logger.setLevel("DEBUG")

    return logger


def attach_file_handler(
    logger: logging.Logger,
    log_path: str,
    log_level: LOG_LEVEL_TYPE = "INFO",
) -> logging.FileHandler:
    """Attach a file handler to a logger."""
    formatter = logging.Formatter("%(asctime)s [%(levelname)s] %(message)s")
    file_handler = logging.FileHandler(log_path)
    file_handler.setLevel(logging.getLevelName(log_level))
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)
    return file_handler


@contextmanager
def file_handler_context(
    logger: logging.Logger,
    log_path: str,
    log_level: LOG_LEVEL_TYPE = "INFO",
):
    """Temporarily attach a file handler to a logger."""
    file_handler = attach_file_handler(logger, log_path, log_level)
    try:
        yield
    finally:
        logger.removeHandler(file_handler)


def set_stream_level(logger: logging.Logger, log_level: LOG_LEVEL_TYPE):
    """Set the levels of stream handlers of a logger."""
    for handler in logger.handlers:
        if isinstance(handler, logging.StreamHandler):
            handler.setLevel(logging.getLevelName(log_level))


@contextmanager
def stream_level_context(logger: logging.Logger, log_level: LOG_LEVEL_TYPE):
    """Temporarily set the levels of stream handlers of a logger."""
    handler_levels = []
    for handler in logger.handlers:
        if isinstance(handler, logging.StreamHandler):
            handler_levels.append((handler, handler.level))
            handler.setLevel(logging.getLevelName(log_level))
    try:
        yield
    finally:
        for handler, orig_level in handler_levels:
            handler.setLevel(orig_level)


@contextmanager
def log_level_context(logger: logging.Logger, log_level: LOG_LEVEL_TYPE):
    """Temporarily set a log level."""
    prev_log_level = logger.level
    logger.setLevel(logging.getLevelName(log_level))
    try:
        yield
    finally:
        logger.setLevel(prev_log_level)
