"""
Functions for importing ecological data.

This module provides functions to import:
- OTU/ASV tables
- Taxonomic classifications
- Sample metadata
- Phylogenetic trees
"""

import pandas as pd
import numpy as np
import biom
import dendropy
from .data_structures import EcologicalData


def read_abundance_table(file_path, format="tabular", transpose=False, **kwargs):
    """
    Read abundance data from file.
    
    Parameters:
        file_path (str): Path to the abundance data file
        format (str): File format:
            - "tabular": Tab-delimited text file
            - "biom": BIOM format file
        transpose (bool): Whether to transpose the table (samples as rows)
        **kwargs: Additional arguments to pass to the reader
    
    Returns:
        pd.DataFrame: Abundance data
    """
    if format == "tabular":
        # Read tab-delimited file
        data = pd.read_csv(file_path, sep="\t", **kwargs)
        # Assume first column is row index
        if 'index_col' not in kwargs:
            data = data.set_index(data.columns[0])
    elif format == "biom":
        # Read BIOM format file
        biom_table = biom.load_table(file_path)
        # Convert to pandas DataFrame
        data = pd.DataFrame(
            biom_table.matrix_data.toarray(),
            index=biom_table.ids('observation'),
            columns=biom_table.ids('sample')
        )
    else:
        raise ValueError(f"Unknown format: {format}")
    
    # Transpose if necessary
    if transpose:
        data = data.T
    
    return data


def read_taxonomy(file_path, format="tabular", **kwargs):
    """
    Read taxonomic classification data from file.
    
    Parameters:
        file_path (str): Path to the taxonomy file
        format (str): File format:
            - "tabular": Tab-delimited text file
        **kwargs: Additional arguments to pass to the reader
    
    Returns:
        pd.DataFrame: Taxonomic classifications
    """
    if format == "tabular":
        # Read tab-delimited file
        data = pd.read_csv(file_path, sep="\t", **kwargs)
        # Assume first column is row index
        if 'index_col' not in kwargs:
            data = data.set_index(data.columns[0])
    else:
        raise ValueError(f"Unknown format: {format}")
    
    return data


def read_metadata(file_path, format="tabular", **kwargs):
    """
    Read sample metadata from file.
    
    Parameters:
        file_path (str): Path to the metadata file
        format (str): File format:
            - "tabular": Tab-delimited text file
        **kwargs: Additional arguments to pass to the reader
    
    Returns:
        pd.DataFrame: Sample metadata
    """
    if format == "tabular":
        # Read tab-delimited file
        data = pd.read_csv(file_path, sep="\t", **kwargs)
        # Assume first column is row index
        if 'index_col' not in kwargs:
            data = data.set_index(data.columns[0])
    else:
        raise ValueError(f"Unknown format: {format}")
    
    return data


def read_tree(file_path, format="newick"):
    """
    Read phylogenetic tree from file.
    
    Parameters:
        file_path (str): Path to the tree file
        format (str): File format:
            - "newick": Newick format
            - "nexus": Nexus format
    
    Returns:
        dendropy.Tree: Phylogenetic tree
    """
    if format == "newick":
        # Read Newick format tree
        tree = dendropy.Tree.get(path=file_path, schema="newick")
    elif format == "nexus":
        # Read Nexus format tree
        tree = dendropy.Tree.get(path=file_path, schema="nexus")
    else:
        raise ValueError(f"Unknown format: {format}")
    
    return tree


def import_data(abundance_file=None, taxonomy_file=None, metadata_file=None, tree_file=None,
                abundance_format="tabular", taxonomy_format="tabular", metadata_format="tabular", tree_format="newick",
                **kwargs):
    """
    Import all data and create an EcologicalData object.
    
    Parameters:
        abundance_file (str): Path to abundance data file
        taxonomy_file (str): Path to taxonomy file
        metadata_file (str): Path to metadata file
        tree_file (str): Path to tree file
        abundance_format (str): Format of abundance file
        taxonomy_format (str): Format of taxonomy file
        metadata_format (str): Format of metadata file
        tree_format (str): Format of tree file
        **kwargs: Additional arguments to pass to readers
    
    Returns:
        EcologicalData: Data object containing all imported data
    """
    # Import abundance data
    abundance = None
    if abundance_file is not None:
        abundance = read_abundance_table(abundance_file, format=abundance_format, **kwargs)
    
    # Import taxonomy
    taxonomy = None
    if taxonomy_file is not None:
        taxonomy = read_taxonomy(taxonomy_file, format=taxonomy_format, **kwargs)
    
    # Import metadata
    metadata = None
    if metadata_file is not None:
        metadata = read_metadata(metadata_file, format=metadata_format, **kwargs)
    
    # Import tree
    tree = None
    if tree_file is not None:
        tree = read_tree(tree_file, format=tree_format)
    
    # Create EcologicalData object
    return EcologicalData(abundance=abundance, taxonomy=taxonomy, metadata=metadata, tree=tree)