# Gene-Disease Network Analysis

This project creates and analyzes graph representations of gene-disease associations retrieved from the DisGeNET database API. The network allows for visualization and analysis of relationships between genes and diseases, with a focus on neurodegenerative disorders.

## Overview

The project includes several components:

- Data retrieval from DisGeNET API
- Construction of gene-disease association graphs
- Multi-gene analysis for deeper insights
- Community detection to identify clusters in the graph
- Visualization of networks and communities

## Getting Started

1. Create a `.env` file with your DisGeNET API key:

```
DISGENET_API_KEY = "your_api_key_here"
```

2. Install dependencies:

```bash
pip install requests networkx matplotlib python-dotenv python-louvain pandas
```

3. Run the scripts in sequence:

```bash
# Collect data for a single gene (APP)
python collectData.py

# Create a graph for a single gene
python create_graph.py

# Create a graph for multiple genes
python multi_gene_graph.py

# Analyze communities in the multi-gene graph
python analyze_communities.py
```

## Key Insights

- Each gene forms its own community with a specific set of associated diseases
- Some diseases (notably Alzheimer's) are linked to multiple genes
- The graph structure reveals patterns in gene-disease relationships

## Example Results

When analyzing five genes (APP, MAPT, PSEN1, SNCA, BACE1) associated with neurodegenerative diseases:

- Alzheimer's Disease is associated with all five genes
- Multiple disease connections exist between specific gene pairs
- Graph community detection identifies gene-specific clusters

For detailed documentation, code explanations, and analysis results, see [DOCUMENTATION.md](./DOCUMENTATION.md).
