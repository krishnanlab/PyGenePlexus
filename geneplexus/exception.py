"""GenePlexus exceptions."""


class CustomNetworkError(Exception):
    """Raised when custom network data files not set up correctly."""


class DownloadError(Exception):
    """Raised when failed to download from URL within max number of retries."""
