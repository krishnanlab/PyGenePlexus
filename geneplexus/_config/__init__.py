"""Basic configurations."""
from . import config
from .logger_util import make_logger

logger = make_logger("INFO")

__all__ = ["config", "logger"]
