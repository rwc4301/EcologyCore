# EcologyCore (Python)

A Python package for ecological data analysis, with a focus on microbial ecology. This is a Python implementation of the original R-based EcologyCore package.

## Features

- Data import and preprocessing for abundance data, taxonomy, metadata, and phylogenetic trees
- Alpha diversity analysis (richness, Shannon, Simpson, etc.)
- Beta diversity analysis (Bray-Curtis, UniFrac, weighted UniFrac)
- Core microbiome identification
- Differential expression analysis
- Compositional regression
- Environmental filtering analysis
- Community assembly processes

## Installation

There are multiple ways to install EcologyCore:

### Option 1: Using the install script (recommended)

```bash
# Clone the repository
git clone https://github.com/username/ecology-core.git
cd ecology-core

# Run the installation script
./install.sh
```

The script provides options to install with conda, pip, or uv based on what's available on your system.

### Option 2: With conda

```bash
# Clone the repository
git clone https://github.com/username/ecology-core.git
cd ecology-core

# Create conda environment
conda env create -f environment.yml

# Activate the environment
conda activate ecology-core
```

### Option 3: With pip

```bash
# Clone the repository
git clone https://github.com/username/ecology-core.git
cd ecology-core

# Create a virtual environment
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate

# Install dependencies
pip install -r requirements.txt

# Install in development mode
pip install -e .
```

## Dependencies

- Python 3.9+
- numpy
- pandas
- matplotlib
- seaborn
- scipy
- scikit-bio
- dendropy
- biom-format
- tqdm

## Getting Started

```python
import numpy as np
import pandas as pd
from ecology_core.data_structures import EcologicalData
from ecology_core.alpha_diversity import calculate_alpha_diversity
from ecology_core.beta_diversity import calculate_beta_diversity, ordinate
from ecology_core.plotting import alpha_diversity_boxplot, ordination_plot

# Create or load your data
abundance = pd.DataFrame(...)  # OTU table (features x samples)
taxonomy = pd.DataFrame(...)   # Taxonomic classifications
metadata = pd.DataFrame(...)   # Sample metadata

# Create EcologicalData object
data = EcologicalData(abundance=abundance, taxonomy=taxonomy, metadata=metadata)

# Preprocess data
filtered_data = data.filter_features(min_prevalence=0.1, min_abundance=10)
rarefied_data = filtered_data.rarefy(depth=10000, seed=42)

# Calculate alpha diversity
alpha_div = calculate_alpha_diversity(rarefied_data, metrics=["observed", "shannon", "simpson"])

# Calculate beta diversity
beta_div = calculate_beta_diversity(rarefied_data, metric="bray")
pcoa = ordinate(beta_div, method="pcoa")

# Create visualizations
alpha_plot = alpha_diversity_boxplot(alpha_div, metric="shannon", group_col="Group")
beta_plot = ordination_plot(pcoa, metadata=metadata, color_col="Group")
```

See the `examples` directory for complete workflows.

## Module Structure

- `data_structures.py`: Core data structures for ecological data
- `import_data.py`: Functions for data import and preprocessing
- `alpha_diversity.py`: Within-sample diversity metrics
- `beta_diversity.py`: Between-sample diversity metrics
- `beta_dispersion.py`: Analysis of community dispersion
- `core_microbiome.py`: Core microbiome identification
- `differential_expression.py`: Differential abundance analysis
- `nst.py`: Normalized Stochasticity Ratio analysis
- `compositional_regression.py`: Regression tools for compositional data
- `environmental_filtering.py`: Analysis of environmental filtering
- `raup_crick.py`: Raup-Crick dissimilarity metrics
- `sncm_fit.py`: Neutral community model fitting
- `plotting.py`: Visualization tools

## Comparison with R Version

This Python implementation aims to provide the same functionality as the R version but with idiomatic Python code and integration with the Python data science ecosystem. Key differences include:

1. **Data Structure**: Uses pandas DataFrames and DendroPy trees instead of phyloseq objects
2. **Visualization**: Uses matplotlib and seaborn instead of ggplot2
3. **Ecosystem Integration**: Leverages scikit-bio for diversity metrics instead of vegan

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.