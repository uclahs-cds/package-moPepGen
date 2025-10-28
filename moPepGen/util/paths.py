"""
Path utilities for moPepGen file operations.

This module provides helper functions for constructing standardized file paths
used throughout the moPepGen pipeline.
"""
from pathlib import Path


def get_peptide_table_path(path: Path) -> Path:
    """
    Get the peptide table file path from the given peptide FASTA path.

    Args:
        path: Path to the peptide FASTA file

    Returns:
        Path to the corresponding peptide table file

    Example:
        >>> get_peptide_table_path(Path('/data/output.fasta'))
        PosixPath('/data/output_peptide_table.txt')
    """
    return path.parent / f"{path.stem}_peptide_table.txt"


def get_peptide_table_path_temp(path: Path) -> Path:
    """
    Get the temporary peptide table path for the given peptide FASTA path.

    Args:
        path: Path to the peptide FASTA file

    Returns:
        Path to the temporary peptide table file

    Example:
        >>> get_peptide_table_path_temp(Path('/data/output.fasta'))
        PosixPath('/data/output_peptide_table.tmp')
    """
    return path.parent / f"{path.stem}_peptide_table.tmp"
