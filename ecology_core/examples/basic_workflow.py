"""
Basic workflow example for EcologyCore.

This script demonstrates a simple analysis workflow using the EcologyCore package:
1. Import synthetic data
2. Preprocess data (filtering, rarefaction)
3. Calculate alpha diversity
4. Calculate beta diversity and ordination
5. Identify core microbiome
6. Create visualizations
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from ecology_core.data_structures import EcologicalData
from ecology_core.alpha_diversity import calculate_alpha_diversity, alpha_diversity_stats
from ecology_core.beta_diversity import calculate_beta_diversity, ordinate
from ecology_core.core_microbiome import identify_core_microbiome
from ecology_core.plotting import alpha_diversity_boxplot, ordination_plot, taxonomic_barplot

# Create a directory for saving plots
os.makedirs("plots", exist_ok=True)

# ============================================================================
# Generate synthetic data for demonstration
# ============================================================================

def generate_synthetic_data(n_samples=30, n_features=100, n_groups=3):
    """Generate synthetic microbiome data for demonstration."""
    # Create sample groups
    group_sizes = [n_samples // n_groups] * n_groups
    group_sizes[-1] += n_samples % n_groups  # Adjust last group if needed
    
    # Create sample IDs and groups
    sample_ids = [f"S{i+1}" for i in range(n_samples)]
    sample_groups = []
    for i, size in enumerate(group_sizes):
        sample_groups.extend([f"Group{i+1}"] * size)
    
    # Create metadata
    metadata = pd.DataFrame({
        'Group': sample_groups,
        'Var1': np.random.normal(0, 1, n_samples),
        'Var2': np.random.uniform(0, 10, n_samples)
    }, index=sample_ids)
    
    # Create feature IDs
    feature_ids = [f"OTU_{i+1}" for i in range(n_features)]
    
    # Generate abundance data with group differences
    abundance = np.zeros((n_features, n_samples))
    
    # Background abundance for all groups
    background = np.random.negative_binomial(5, 0.5, (n_features, n_samples))
    
    # Group-specific patterns
    start_idx = 0
    for i, size in enumerate(group_sizes):
        end_idx = start_idx + size
        
        # Features that are enriched in this group
        enriched_features = np.random.choice(range(n_features), size=15, replace=False)
        
        # Add enrichment
        abundance[enriched_features, start_idx:end_idx] += np.random.negative_binomial(20, 0.3, (15, size))
        
        start_idx = end_idx
    
    # Add background and convert to DataFrame
    abundance = pd.DataFrame(abundance + background, index=feature_ids, columns=sample_ids)
    
    # Create taxonomy data
    taxa = ['Firmicutes', 'Bacteroidetes', 'Proteobacteria', 'Actinobacteria', 'Verrucomicrobia']
    phylum = np.random.choice(taxa, n_features)
    
    genera = ['Genus_' + str(i) for i in range(1, 21)]
    genus = np.random.choice(genera, n_features)
    
    taxonomy = pd.DataFrame({
        'Kingdom': 'Bacteria',
        'Phylum': phylum,
        'Class': [f"Class_{i%10+1}" for i in range(n_features)],
        'Order': [f"Order_{i%15+1}" for i in range(n_features)],
        'Family': [f"Family_{i%25+1}" for i in range(n_features)],
        'Genus': genus,
        'Species': [f"Species_{i+1}" for i in range(n_features)]
    }, index=feature_ids)
    
    # Create EcologicalData object
    data = EcologicalData(abundance=abundance, taxonomy=taxonomy, metadata=metadata)
    
    return data

# Generate synthetic data
print("Generating synthetic data...")
data = generate_synthetic_data()
print(f"Generated data with {data.abundance.shape[0]} features and {data.abundance.shape[1]} samples")

# ============================================================================
# Data preprocessing
# ============================================================================

# Filter low-abundance features
print("Filtering low-abundance features...")
filtered_data = data.filter_features(min_prevalence=0.1, min_abundance=10)
print(f"Retained {filtered_data.abundance.shape[0]} features after filtering")

# Rarefy to even sequencing depth
print("Rarefying data...")
rarefied_data = filtered_data.rarefy(seed=42)
print(f"Rarefied to depth of {rarefied_data.abundance.sum(axis=0).iloc[0]} reads per sample")

# ============================================================================
# Alpha diversity analysis
# ============================================================================

# Calculate alpha diversity metrics
print("Calculating alpha diversity...")
alpha_div = calculate_alpha_diversity(rarefied_data, metrics="all")
print(f"Calculated alpha diversity metrics: {', '.join(alpha_div.columns[:5])}")

# Test for differences between groups
print("Testing for differences in Shannon diversity between groups...")
shannon_stats = alpha_diversity_stats(alpha_div, group_col="Group", metric="shannon", test="anova")
print(f"ANOVA results: F={shannon_stats['statistic']:.2f}, p={shannon_stats['p_value']:.4f}")

# Create alpha diversity boxplot
print("Creating alpha diversity boxplot...")
fig = alpha_diversity_boxplot(alpha_div, metric="shannon", group_col="Group")
fig.savefig("plots/alpha_diversity_shannon.png", dpi=300)
plt.close(fig)

# ============================================================================
# Beta diversity analysis
# ============================================================================

# Calculate beta diversity
print("Calculating beta diversity (Bray-Curtis)...")
beta_div = calculate_beta_diversity(rarefied_data, metric="bray", binary=False)

# Perform ordination
print("Performing PCoA ordination...")
ordination_results = ordinate(beta_div, method="pcoa", n_components=2)

# Create ordination plot
print("Creating ordination plot...")
fig = ordination_plot(
    ordination_results, 
    metadata=rarefied_data.metadata, 
    color_col="Group",
    title="PCoA of Bray-Curtis Distances"
)
fig.savefig("plots/pcoa_bray_curtis.png", dpi=300)
plt.close(fig)

# ============================================================================
# Core microbiome analysis
# ============================================================================

# Identify core microbiome for each group
print("Identifying core microbiome...")
core_microbiome = identify_core_microbiome(
    rarefied_data,
    prevalence_threshold=0.7,
    abundance_threshold=0.001,
    group_col="Group"
)

for group, core_features in core_microbiome.items():
    print(f"Group {group} core microbiome: {len(core_features)} features")

# ============================================================================
# Taxonomic composition
# ============================================================================

# Create taxonomic barplot
print("Creating taxonomic composition barplot...")
fig = taxonomic_barplot(
    rarefied_data,
    tax_level='Phylum',
    top_n=5,
    figsize=(10, 6)
)
fig.savefig("plots/taxonomic_composition.png", dpi=300)
plt.close(fig)

print("Analysis completed! Check the 'plots' directory for output visualizations.")