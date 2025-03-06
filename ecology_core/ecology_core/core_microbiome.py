"""
Core microbiome analysis functions.

This module provides functions to identify core microbiome members
based on prevalence and abundance thresholds.
"""

import numpy as np
import pandas as pd
from .data_structures import EcologicalData


def identify_core_microbiome(data, prevalence_threshold=0.8, abundance_threshold=0.0,
                           group_col=None, groups=None):
    """
    Identify core microbiome members based on prevalence and abundance thresholds.
    
    Parameters:
        data (EcologicalData): Ecological data object
        prevalence_threshold (float): Minimum proportion of samples a feature must be present in
        abundance_threshold (float): Minimum relative abundance a feature must have
        group_col (str): Metadata column to use for grouping samples
        groups (list): Specific groups to compute core microbiome for
    
    Returns:
        dict: Core microbiome results for each group, with features and their prevalence/abundance
    """
    # Ensure data is an EcologicalData object
    if not isinstance(data, EcologicalData):
        raise TypeError("data must be an EcologicalData object")
    
    # Calculate relative abundance
    rel_abundance = data.transform(method="relative").abundance
    
    # Initialize results
    core_results = {}
    
    # If no group_col, calculate for all samples
    if group_col is None:
        # Calculate prevalence (proportion of samples with feature present)
        prevalence = (data.abundance > 0).mean(axis=1)
        
        # Calculate mean relative abundance
        mean_abundance = rel_abundance.mean(axis=1)
        
        # Identify core features
        core_features = (prevalence >= prevalence_threshold) & (mean_abundance >= abundance_threshold)
        
        # Create results dataframe
        core_df = pd.DataFrame({
            'prevalence': prevalence[core_features],
            'mean_abundance': mean_abundance[core_features]
        }).sort_values('prevalence', ascending=False)
        
        core_results['all'] = core_df
    
    # If group_col is provided, calculate core microbiome for each group
    else:
        # Check if metadata is available
        if data.metadata is None:
            raise ValueError("Metadata is required for group-specific core microbiome")
        
        # Check if group_col exists in metadata
        if group_col not in data.metadata.columns:
            raise ValueError(f"Group column '{group_col}' not found in metadata")
        
        # Get groups to analyze
        if groups is None:
            groups = data.metadata[group_col].unique()
        
        # Calculate core microbiome for each group
        for group in groups:
            # Get samples in this group
            group_samples = data.metadata[data.metadata[group_col] == group].index
            
            # Skip if no samples
            if len(group_samples) == 0:
                continue
            
            # Filter abundance data to group samples
            group_abundance = data.abundance[group_samples]
            group_rel_abundance = rel_abundance[group_samples]
            
            # Calculate prevalence and mean abundance
            prevalence = (group_abundance > 0).mean(axis=1)
            mean_abundance = group_rel_abundance.mean(axis=1)
            
            # Identify core features
            core_features = (prevalence >= prevalence_threshold) & (mean_abundance >= abundance_threshold)
            
            # Create results dataframe
            core_df = pd.DataFrame({
                'prevalence': prevalence[core_features],
                'mean_abundance': mean_abundance[core_features]
            }).sort_values('prevalence', ascending=False)
            
            core_results[group] = core_df
    
    return core_results


def core_microbiome_venn(data, group_col, groups=None, prevalence_threshold=0.8, 
                        abundance_threshold=0.0):
    """
    Calculate overlapping core microbiome features between groups for Venn diagrams.
    
    Parameters:
        data (EcologicalData): Ecological data object
        group_col (str): Metadata column to use for grouping samples
        groups (list): Specific groups to include (up to 5)
        prevalence_threshold (float): Minimum proportion of samples a feature must be present in
        abundance_threshold (float): Minimum relative abundance a feature must have
    
    Returns:
        dict: Sets of core features for each group and their intersections
    """
    # Ensure data is an EcologicalData object
    if not isinstance(data, EcologicalData):
        raise TypeError("data must be an EcologicalData object")
    
    # Check if metadata is available
    if data.metadata is None:
        raise ValueError("Metadata is required for group-specific core microbiome")
    
    # Check if group_col exists in metadata
    if group_col not in data.metadata.columns:
        raise ValueError(f"Group column '{group_col}' not found in metadata")
    
    # Get groups to analyze
    if groups is None:
        groups = data.metadata[group_col].unique()
    
    # Limit to 5 groups for Venn diagram
    if len(groups) > 5:
        raise ValueError("Maximum 5 groups supported for Venn diagram")
    
    # Get core microbiome for each group
    core_results = identify_core_microbiome(
        data, 
        prevalence_threshold=prevalence_threshold,
        abundance_threshold=abundance_threshold,
        group_col=group_col,
        groups=groups
    )
    
    # Extract core feature sets
    group_features = {group: set(df.index) for group, df in core_results.items()}
    
    # Calculate set intersections
    intersections = {}
    
    # Single group sets
    for i, group in enumerate(groups):
        if group in group_features:
            set_name = str(group)
            intersections[set_name] = group_features[group]
    
    # Pair intersections
    for i, group1 in enumerate(groups):
        for j, group2 in enumerate(groups[i+1:], i+1):
            if group1 in group_features and group2 in group_features:
                set_name = f"{group1} & {group2}"
                intersections[set_name] = group_features[group1].intersection(group_features[group2])
    
    # Triple intersections
    for i, group1 in enumerate(groups):
        for j, group2 in enumerate(groups[i+1:], i+1):
            for k, group3 in enumerate(groups[j+1:], j+1):
                if all(g in group_features for g in [group1, group2, group3]):
                    set_name = f"{group1} & {group2} & {group3}"
                    intersections[set_name] = group_features[group1].intersection(
                        group_features[group2], group_features[group3]
                    )
    
    # Quadruple intersections
    if len(groups) >= 4:
        for i, group1 in enumerate(groups):
            for j, group2 in enumerate(groups[i+1:], i+1):
                for k, group3 in enumerate(groups[j+1:], j+1):
                    for l, group4 in enumerate(groups[k+1:], k+1):
                        if all(g in group_features for g in [group1, group2, group3, group4]):
                            set_name = f"{group1} & {group2} & {group3} & {group4}"
                            intersections[set_name] = group_features[group1].intersection(
                                group_features[group2], group_features[group3], group_features[group4]
                            )
    
    # Quintuple intersection
    if len(groups) == 5:
        if all(g in group_features for g in groups):
            set_name = " & ".join(groups)
            intersections[set_name] = set.intersection(*[group_features[g] for g in groups])
    
    return intersections


def get_core_taxonomy(data, core_features, tax_level='Phylum'):
    """
    Summarize taxonomy of core microbiome features.
    
    Parameters:
        data (EcologicalData): Ecological data object
        core_features (list or set): Core microbiome features
        tax_level (str): Taxonomic level to summarize
    
    Returns:
        pd.DataFrame: Taxonomic summary of core features
    """
    # Ensure data is an EcologicalData object
    if not isinstance(data, EcologicalData):
        raise TypeError("data must be an EcologicalData object")
    
    # Check if taxonomy data is available
    if data.taxonomy is None:
        raise ValueError("Taxonomy data is required for core taxonomy summary")
    
    # Check if tax_level exists in taxonomy
    if tax_level not in data.taxonomy.columns:
        raise ValueError(f"Taxonomic level '{tax_level}' not found in taxonomy data")
    
    # Create a set of core features
    core_set = set(core_features)
    
    # Extract taxonomy for core features
    core_tax = data.taxonomy.loc[data.taxonomy.index.intersection(core_set)]
    
    # Summarize taxonomy
    tax_counts = core_tax[tax_level].value_counts()
    tax_percent = (tax_counts / len(core_tax)) * 100
    
    # Create summary dataframe
    tax_summary = pd.DataFrame({
        'count': tax_counts,
        'percent': tax_percent
    })
    
    return tax_summary