"""GenePlexus exceptions."""


class CustomDataError(Exception):
    """Raised when custom network or gsc data files not set up correctly."""


class DownloadError(Exception):
    """Raised when failed to download from URL within max number of retries."""
