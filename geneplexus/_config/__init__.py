"""Basic configurations."""
from . import config
from .logger_util import make_logger

logger = make_logger(log_level="INFO")

__all__ = ["config", "logger"]
