from pkg_resources import get_distribution

try:
    __version__ = get_distribution("viridian").version
except:
    __version__ = "local"


__all__ = [
    "amplicon_overlapper",
    "amplicons",
    "assemble",
    "racon",
    "tasks",
    "utils",
]

from viridian import *
