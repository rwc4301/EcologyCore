"""
Core data structures for representing ecological data.

This module provides classes for storing and manipulating:
- Abundance data (OTU/ASV tables)
- Taxonomic information
- Sample metadata
- Phylogenetic trees
"""

import pandas as pd
import numpy as np
import dendropy
from scipy import sparse


class EcologicalData:
    """
    Main container class for ecological data, similar to phyloseq in R.
    
    Attributes:
        abundance (pd.DataFrame): OTU/ASV abundance table (features x samples)
        taxonomy (pd.DataFrame): Taxonomic classification for each feature
        metadata (pd.DataFrame): Sample metadata
        tree (dendropy.Tree): Phylogenetic tree
    """
    
    def __init__(self, abundance=None, taxonomy=None, metadata=None, tree=None):
        """
        Initialize an EcologicalData object.
        
        Parameters:
            abundance (pd.DataFrame): OTU/ASV abundance table (features x samples)
            taxonomy (pd.DataFrame): Taxonomic classification for each feature
            metadata (pd.DataFrame): Sample metadata
            tree (dendropy.Tree): Phylogenetic tree
        """
        self.abundance = abundance
        self.taxonomy = taxonomy
        self.metadata = metadata
        self.tree = tree
        self._validate()
    
    def _validate(self):
        """Validate the data structures for consistency."""
        if self.abundance is not None:
            if self.taxonomy is not None:
                # Check that all OTUs in abundance have taxonomy
                missing = set(self.abundance.index) - set(self.taxonomy.index)
                if missing:
                    print(f"Warning: {len(missing)} OTUs in abundance data don't have taxonomy")
            
            if self.metadata is not None:
                # Check that all samples in abundance have metadata
                missing = set(self.abundance.columns) - set(self.metadata.index)
                if missing:
                    print(f"Warning: {len(missing)} samples in abundance data don't have metadata")
    
    def filter_features(self, min_prevalence=0, min_abundance=0):
        """
        Filter features based on prevalence and abundance.
        
        Parameters:
            min_prevalence (float): Minimum proportion of samples a feature must be present in
            min_abundance (float): Minimum total abundance a feature must have
        
        Returns:
            EcologicalData: Filtered data object
        """
        if self.abundance is None:
            return self
        
        # Calculate prevalence
        prevalence = (self.abundance > 0).mean(axis=1)
        # Calculate total abundance
        total_abundance = self.abundance.sum(axis=1)
        
        # Filter features
        keep = (prevalence >= min_prevalence) & (total_abundance >= min_abundance)
        filtered_abundance = self.abundance.loc[keep]
        
        # Filter taxonomy if present
        filtered_taxonomy = None
        if self.taxonomy is not None:
            filtered_taxonomy = self.taxonomy.loc[filtered_taxonomy.index.intersection(filtered_abundance.index)]
        
        # Return new object
        return EcologicalData(
            abundance=filtered_abundance,
            taxonomy=filtered_taxonomy,
            metadata=self.metadata,
            tree=self.tree  # Tree filtering not implemented
        )
    
    def rarefy(self, depth=None, seed=None):
        """
        Rarefy abundance data to a specified sequencing depth.
        
        Parameters:
            depth (int): Sequencing depth for rarefaction
            seed (int): Random seed for reproducibility
        
        Returns:
            EcologicalData: Rarefied data object
        """
        if self.abundance is None:
            return self
        
        if depth is None:
            # Default to minimum sample sum
            depth = self.abundance.sum(axis=0).min()
        
        np.random.seed(seed)
        
        # Function to rarefy a single sample
        def rarefy_sample(x):
            if sum(x) <= depth:
                return x
            prob = x / sum(x)
            return np.random.multinomial(depth, prob)
        
        # Apply rarefaction to each sample
        rarefied = self.abundance.apply(rarefy_sample, axis=0)
        
        # Return new object
        return EcologicalData(
            abundance=rarefied,
            taxonomy=self.taxonomy,
            metadata=self.metadata,
            tree=self.tree
        )
    
    def transform(self, method="relative"):
        """
        Transform abundance data.
        
        Parameters:
            method (str): Transformation method:
                - "relative": Convert to relative abundance
                - "log": Log transformation (log(x+1))
                - "clr": Centered log-ratio transformation
        
        Returns:
            EcologicalData: Transformed data object
        """
        if self.abundance is None:
            return self
        
        if method == "relative":
            # Convert to relative abundance
            transformed = self.abundance.div(self.abundance.sum(axis=0), axis=1)
        elif method == "log":
            # Log transformation
            transformed = np.log1p(self.abundance)
        elif method == "clr":
            # Centered log-ratio transformation
            # Add pseudocount to avoid log(0)
            transformed = self.abundance + 1
            # Calculate geometric mean of each sample
            geometric_mean = transformed.apply(lambda x: np.exp(np.mean(np.log(x))), axis=0)
            # Apply CLR transformation
            transformed = transformed.apply(lambda x: np.log(x / np.exp(np.mean(np.log(x)))), axis=0)
        else:
            raise ValueError(f"Unknown transformation method: {method}")
        
        # Return new object
        return EcologicalData(
            abundance=transformed,
            taxonomy=self.taxonomy,
            metadata=self.metadata,
            tree=self.tree
        )