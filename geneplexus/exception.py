"""GenePlexus exceptions."""


class FlyMonarchError(Exception):
    """Raised becasue no Monarch annotations for Fly."""

    
class MondoError(Exception):
    """Raised becasue no Mondo only has annotations for Human."""


class ZebrafishBioGRIDError(Exception):
    """Raised becuse BioGRID doesn't have a good Zebrafish network."""


class NoPositivesError(Exception):
    """Raised when there are no positive genes in the network."""


class DownloadError(Exception):
    """Raised when failed to download from URL within max number of retries."""
