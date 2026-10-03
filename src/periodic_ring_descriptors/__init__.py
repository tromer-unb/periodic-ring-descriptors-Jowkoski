"""Periodic ring descriptors for 2D and 3D network materials."""

from .descriptor_2d import compute_descriptor_2d, joukowsky_metrics_v2
from .descriptor_3d import (
    DESCRIPTOR_VERSION as DESCRIPTOR_3D_VERSION,
    compute_descriptor_3d,
    compute_descriptor_3d_structure,
)
from .features import structure_feature_vector

__all__ = [
    "compute_descriptor_2d",
    "compute_descriptor_3d",
    "compute_descriptor_3d_structure",
    "joukowsky_metrics_v2",
    "structure_feature_vector",
    "DESCRIPTOR_3D_VERSION",
]

__version__ = "0.1.0"