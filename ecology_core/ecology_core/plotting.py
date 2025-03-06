"""
Visualization functions for ecological data.

This module provides functions to create common plots for ecological data analysis:
- Alpha diversity boxplots
- Beta diversity ordination plots
- Taxonomic composition barplots
- Heatmaps
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from .data_structures import EcologicalData


def alpha_diversity_boxplot(alpha_div, metric, group_col, palette=None, title=None, figsize=(10, 6)):
    """
    Create a boxplot of alpha diversity metrics grouped by a metadata column.
    
    Parameters:
        alpha_div (pd.DataFrame): Alpha diversity metrics from calculate_alpha_diversity
        metric (str): Column name of the metric to plot
        group_col (str): Column name for grouping samples
        palette (str or list): Color palette to use
        title (str): Plot title
        figsize (tuple): Figure size (width, height)
    
    Returns:
        matplotlib.figure.Figure: The figure object
    """
    # Create figure
    fig, ax = plt.subplots(figsize=figsize)
    
    # Create boxplot
    sns.boxplot(
        data=alpha_div,
        x=group_col,
        y=metric,
        palette=palette,
        ax=ax
    )
    
    # Add individual points
    sns.stripplot(
        data=alpha_div,
        x=group_col,
        y=metric,
        color='black',
        size=3,
        alpha=0.5,
        jitter=True,
        ax=ax
    )
    
    # Set title and labels
    if title:
        ax.set_title(title)
    ax.set_xlabel(group_col)
    ax.set_ylabel(metric)
    
    # Rotate x-axis labels if needed
    plt.xticks(rotation=45, ha='right')
    
    plt.tight_layout()
    return fig


def ordination_plot(ordination, metadata=None, color_col=None, shape_col=None, 
                   title=None, figsize=(10, 8)):
    """
    Create an ordination plot (e.g., PCoA, NMDS) with optional metadata overlays.
    
    Parameters:
        ordination (pd.DataFrame): Ordination results from ordinate function
        metadata (pd.DataFrame): Sample metadata
        color_col (str): Metadata column to use for point colors
        shape_col (str): Metadata column to use for point shapes
        title (str): Plot title
        figsize (tuple): Figure size (width, height)
    
    Returns:
        matplotlib.figure.Figure: The figure object
    """
    # Create figure
    fig, ax = plt.subplots(figsize=figsize)
    
    # Prepare plot data
    plot_data = ordination.copy()
    
    # Add metadata if provided
    if metadata is not None and (color_col is not None or shape_col is not None):
        plot_data = pd.merge(
            plot_data, 
            metadata, 
            left_index=True, 
            right_index=True, 
            how="left"
        )
    
    # Extract axis labels
    x_col = plot_data.columns[0]
    y_col = plot_data.columns[1]
    
    # Create scatter plot
    if color_col and shape_col:
        # Plot with both color and shape
        for name, group in plot_data.groupby([color_col, shape_col]):
            ax.scatter(
                group[x_col], 
                group[y_col], 
                label=f"{name[0]}, {name[1]}",
                alpha=0.8
            )
    elif color_col:
        # Plot with color only
        for name, group in plot_data.groupby(color_col):
            ax.scatter(
                group[x_col], 
                group[y_col], 
                label=name,
                alpha=0.8
            )
    elif shape_col:
        # Plot with shape only
        for name, group in plot_data.groupby(shape_col):
            ax.scatter(
                group[x_col], 
                group[y_col], 
                label=name,
                alpha=0.8
            )
    else:
        # Simple scatter plot without grouping
        ax.scatter(plot_data[x_col], plot_data[y_col], alpha=0.8)
    
    # Add axis labels and title
    ax.set_xlabel(x_col)
    ax.set_ylabel(y_col)
    if title:
        ax.set_title(title)
    
    # Add legend if needed
    if color_col or shape_col:
        ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    
    # Add grid
    ax.grid(True, linestyle='--', alpha=0.7)
    
    # Add horizontal and vertical lines at origin
    ax.axhline(y=0, color='k', linestyle='-', alpha=0.3)
    ax.axvline(x=0, color='k', linestyle='-', alpha=0.3)
    
    plt.tight_layout()
    return fig


def taxonomic_barplot(data, tax_level='Phylum', top_n=10, sample_groups=None, 
                    palette=None, figsize=(12, 8)):
    """
    Create a stacked barplot of taxonomic composition.
    
    Parameters:
        data (EcologicalData): Ecological data object
        tax_level (str): Taxonomic level to plot
        top_n (int): Number of top taxa to display, remaining are grouped as 'Other'
        sample_groups (dict): Dictionary mapping sample IDs to groups for organizing
        palette (str or list): Color palette to use
        figsize (tuple): Figure size (width, height)
    
    Returns:
        matplotlib.figure.Figure: The figure object
    """
    # Ensure data is an EcologicalData object
    if not isinstance(data, EcologicalData):
        raise TypeError("data must be an EcologicalData object")
    
    # Check if taxonomy data is available
    if data.taxonomy is None:
        raise ValueError("Taxonomy data is required for taxonomic barplot")
    
    # Transform abundance to relative abundance
    rel_abundance = data.transform(method="relative").abundance
    
    # Group by the specified taxonomic level
    if tax_level not in data.taxonomy.columns:
        raise ValueError(f"Taxonomic level '{tax_level}' not found in taxonomy data")
    
    taxa_groups = {}
    for feature in rel_abundance.index:
        if feature in data.taxonomy.index:
            taxon = data.taxonomy.loc[feature, tax_level]
            if pd.isna(taxon):
                taxon = "Unclassified"
            if taxon not in taxa_groups:
                taxa_groups[taxon] = []
            taxa_groups[taxon].append(feature)
        else:
            # Handle features without taxonomy
            if "Unclassified" not in taxa_groups:
                taxa_groups["Unclassified"] = []
            taxa_groups["Unclassified"].append(feature)
    
    # Sum abundances by taxonomic group
    taxa_abundance = pd.DataFrame(index=taxa_groups.keys(), columns=rel_abundance.columns)
    for taxon, features in taxa_groups.items():
        taxa_abundance.loc[taxon] = rel_abundance.loc[features].sum(axis=0)
    
    # Sort taxa by total abundance and keep top N
    taxa_totals = taxa_abundance.sum(axis=1).sort_values(ascending=False)
    top_taxa = taxa_totals.index[:top_n].tolist()
    
    # Group the rest as 'Other'
    if len(taxa_totals) > top_n:
        other_taxa = taxa_totals.index[top_n:].tolist()
        taxa_abundance.loc['Other'] = taxa_abundance.loc[other_taxa].sum(axis=0)
        taxa_abundance = taxa_abundance.loc[top_taxa + ['Other']]
    else:
        taxa_abundance = taxa_abundance.loc[top_taxa]
    
    # Create figure
    fig, ax = plt.subplots(figsize=figsize)
    
    # Reorder samples if sample_groups is provided
    if sample_groups is not None:
        # Create a new column for sample groups
        sample_order = []
        for group, samples in sample_groups.items():
            for sample in samples:
                if sample in taxa_abundance.columns:
                    sample_order.append(sample)
        
        # Filter to keep only samples in the data
        sample_order = [s for s in sample_order if s in taxa_abundance.columns]
        
        # Add any samples not in sample_groups
        missing_samples = [s for s in taxa_abundance.columns if s not in sample_order]
        sample_order.extend(missing_samples)
        
        # Reorder columns
        taxa_abundance = taxa_abundance[sample_order]
    
    # Plot stacked bar chart
    taxa_abundance.T.plot(
        kind='bar', 
        stacked=True, 
        ax=ax, 
        colormap=palette if palette else None
    )
    
    # Set labels and title
    ax.set_xlabel('Sample')
    ax.set_ylabel('Relative Abundance')
    ax.set_title(f'Taxonomic Composition at {tax_level} Level')
    
    # Rotate x-axis labels for readability
    plt.xticks(rotation=45, ha='right')
    
    # Adjust legend
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(
        handles, labels,
        bbox_to_anchor=(1.05, 1), 
        loc='upper left', 
        title=tax_level
    )
    
    plt.tight_layout()
    return fig


def heatmap(data, metadata=None, sample_cluster=True, feature_cluster=True,
           transform='log', top_n=50, figsize=(12, 10)):
    """
    Create a heatmap of abundance data.
    
    Parameters:
        data (EcologicalData): Ecological data object
        metadata (str or list): Metadata columns to include as annotations
        sample_cluster (bool): Whether to cluster samples
        feature_cluster (bool): Whether to cluster features
        transform (str): Transformation to apply to abundance data:
            - "log": Log transformation
            - "clr": Centered log-ratio
            - "relative": Relative abundance
            - None: No transformation
        top_n (int): Number of top abundant features to include
        figsize (tuple): Figure size (width, height)
    
    Returns:
        matplotlib.figure.Figure: The figure object
    """
    # Ensure data is an EcologicalData object
    if not isinstance(data, EcologicalData):
        raise TypeError("data must be an EcologicalData object")
    
    # Transform data if specified
    if transform:
        transformed_data = data.transform(method=transform).abundance
    else:
        transformed_data = data.abundance.copy()
    
    # Select top features by mean abundance
    if top_n and top_n < transformed_data.shape[0]:
        feature_means = transformed_data.mean(axis=1).sort_values(ascending=False)
        top_features = feature_means.index[:top_n]
        transformed_data = transformed_data.loc[top_features]
    
    # Create figure
    fig, ax = plt.subplots(figsize=figsize)
    
    # Create sample annotations if metadata is provided
    row_colors = None
    if metadata is not None and data.metadata is not None:
        if isinstance(metadata, str):
            metadata = [metadata]
        
        # Extract specified metadata columns
        sample_annotations = data.metadata[metadata].copy()
        
        # Create a color mapping for categorical variables
        cmap_list = ['Set1', 'Set2', 'Set3', 'tab10', 'tab20']
        row_colors = []
        
        for i, col in enumerate(sample_annotations.columns):
            # Use a different colormap for each metadata column
            cmap_name = cmap_list[i % len(cmap_list)]
            cmap = plt.cm.get_cmap(cmap_name)
            
            # Map categorical values to colors
            if sample_annotations[col].dtype == 'object':
                unique_values = sample_annotations[col].unique()
                color_dict = {val: cmap(i / len(unique_values)) 
                            for i, val in enumerate(unique_values)}
                row_colors.append(sample_annotations[col].map(color_dict))
            else:
                # For numerical columns, use a colormap directly
                row_colors.append(sample_annotations[col])
    
    # Create the clustered heatmap
    g = sns.clustermap(
        transformed_data,
        cmap='viridis',
        row_cluster=feature_cluster,
        col_cluster=sample_cluster,
        row_colors=row_colors,
        figsize=figsize,
        xticklabels=True,
        yticklabels=False
    )
    
    # Set title
    if transform:
        plt.suptitle(f'{transform.upper()} Transformed Abundance Heatmap', y=1.02)
    else:
        plt.suptitle('Abundance Heatmap', y=1.02)
    
    # Add colorbar label
    g.ax_heatmap.collections[0].colorbar.set_label(
        'Abundance' if transform is None else f'{transform.upper()} Transformed Abundance'
    )
    
    plt.tight_layout()
    return g.fig