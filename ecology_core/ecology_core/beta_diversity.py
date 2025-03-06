"""
Beta diversity metrics for ecological data.

This module provides functions to calculate various beta diversity metrics,
which measure the diversity between different samples.
"""

import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist, squareform
from skbio.diversity import beta
from skbio.stats.ordination import pcoa
from .data_structures import EcologicalData


def calculate_beta_diversity(data, metric="bray", binary=False, rarefy_depth=None):
    """
    Calculate beta diversity distances between samples.
    
    Parameters:
        data (EcologicalData): Ecological data object
        metric (str): Distance metric:
            - "bray": Bray-Curtis dissimilarity
            - "jaccard": Jaccard distance
            - "unifrac": UniFrac distance (requires tree)
            - "wunifrac": Weighted UniFrac distance (requires tree)
            - "euclidean": Euclidean distance
            - "manhattan": Manhattan distance
        binary (bool): Convert abundance to presence/absence
        rarefy_depth (int): Depth for rarefaction before calculating diversity
            Set to None to use raw counts
    
    Returns:
        pd.DataFrame: Distance matrix between samples
    """
    # Ensure data is an EcologicalData object
    if not isinstance(data, EcologicalData):
        raise TypeError("data must be an EcologicalData object")
    
    # Rarefy if specified
    if rarefy_depth is not None:
        data = data.rarefy(depth=rarefy_depth)
    
    # Get abundance data
    abundance = data.abundance.copy()
    
    # Convert to binary if specified
    if binary:
        abundance = (abundance > 0).astype(int)
    
    # Calculate distances based on the specified metric
    if metric == "bray":
        # Bray-Curtis dissimilarity
        distances = beta.bray_curtis(abundance.T)
    elif metric == "jaccard":
        # Jaccard distance
        distances = beta.jaccard(abundance.T)
    elif metric in ["unifrac", "wunifrac"]:
        # UniFrac distance (requires tree)
        if data.tree is None:
            raise ValueError(f"{metric} requires a phylogenetic tree")
        
        # TODO: Implement UniFrac using scikit-bio
        # This is a placeholder - UniFrac implementation would require
        # additional code to map OTUs to tree tips and compute distances
        raise NotImplementedError(f"{metric} is not yet implemented")
    else:
        # Use scipy's pdist for other metrics
        try:
            distances = pdist(abundance.T, metric=metric)
            distances = squareform(distances)
        except ValueError:
            raise ValueError(f"Unknown metric: {metric}")
    
    # Convert to DataFrame
    dist_df = pd.DataFrame(
        distances,
        index=abundance.columns,
        columns=abundance.columns
    )
    
    return dist_df


def ordinate(distances, method="pcoa", n_components=2):
    """
    Perform ordination on beta diversity distances.
    
    Parameters:
        distances (pd.DataFrame): Distance matrix from calculate_beta_diversity
        method (str): Ordination method:
            - "pcoa": Principal Coordinates Analysis
            - "nmds": Non-metric Multidimensional Scaling (not implemented yet)
        n_components (int): Number of components to return
    
    Returns:
        pd.DataFrame: Ordination results with samples as rows and components as columns
    """
    if method == "pcoa":
        # Perform PCoA
        pcoa_results = pcoa(distances)
        
        # Extract the coordinates
        coords = pcoa_results.samples.iloc[:, :n_components]
        
        # Add variance explained as column names
        explained_var = pcoa_results.proportion_explained.iloc[:n_components]
        columns = [f"PCo{i+1} ({var:.2%})" for i, var in enumerate(explained_var)]
        coords.columns = columns
        
        return coords
    
    elif method == "nmds":
        # NMDS is not yet implemented
        raise NotImplementedError("NMDS is not yet implemented")
    
    else:
        raise ValueError(f"Unknown ordination method: {method}")


def beta_permanova(distances, metadata, formula, permutations=999):
    """
    Perform PERMANOVA test on beta diversity distances.
    
    Parameters:
        distances (pd.DataFrame): Distance matrix from calculate_beta_diversity
        metadata (pd.DataFrame): Sample metadata
        formula (str): Formula for the test (e.g., "~ group")
        permutations (int): Number of permutations
    
    Returns:
        pd.DataFrame: PERMANOVA results
    """
    # This requires the implementation of PERMANOVA
    # which is complex and would normally use a package like vegan in R
    # In Python, scikit-bio has a permanova implementation but requires
    # additional dependencies
    
    # This is a placeholder for the actual implementation
    raise NotImplementedError("PERMANOVA is not yet implemented")
    
    # Rough sketch of implementation would be:
    # 1. Create a patsy design matrix from formula and metadata
    # 2. Calculate pseudo-F statistic between groups
    # 3. Permute sample labels and calculate null distribution
    # 4. Calculate p-values from permutations
    
    # Example pseudocode:
    # import patsy
    # from skbio.stats.distance import permanova
    #
    # # Create design matrix
    # design = patsy.dmatrix(formula, metadata, return_type='dataframe')
    #
    # # Run PERMANOVA (if skbio had this implemented)
    # result = permanova(distances.values, design)
    #
    # return result