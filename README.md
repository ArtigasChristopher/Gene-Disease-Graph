# Gene-Disease Network Analysis with Advanced Validation

This project creates and analyzes graph representations of gene-disease associations retrieved from the DisGeNET database API. The network allows for visualization and analysis of relationships between genes and diseases, with a focus on neurodegenerative disorders, including comprehensive validation metrics and statistical analysis.

## Overview

The project includes several components:

- Data retrieval from DisGeNET API
- Construction of gene-disease association graphs
- Multi-gene analysis for deeper insights
- Community detection to identify clusters in the graph
- **Advanced validation with precision metrics (NEW)**
- **Statistical significance testing (NEW)**
- **Network topology validation (NEW)**
- Visualization of networks and communities

## Getting Started

1. Create a `.env` file with your DisGeNET API key:

```
DISGENET_API_KEY = "your_api_key_here"
```

2. Install dependencies:

```bash
pip install requests networkx matplotlib python-dotenv python-louvain pandas scikit-learn seaborn scipy
```

3. Run the scripts in sequence:

```bash
# Collect data for a single gene (APP)
python3 collectData.py

# Create a graph for a single gene
python3 create_graph.py

# Create a graph for multiple genes
python3 multi_gene_graph.py

# Analyze communities in the multi-gene graph
python3 analyze_communities.py

# Run advanced validation and precision analysis (NEW)
python3 advanced_validation.py

# Generate comprehensive summary report (NEW)
python3 generate_summary.py
```

## Validation and Precision Metrics

Our analysis includes comprehensive validation to ensure data quality and reliability:

### Quality Metrics

- **High-confidence associations**: 30.3% (73/241) with strong evidence support
- **Evidence support**: Mean of 173.9 PMID citations per association
- **Score distribution**: Mean 0.598 ± 0.208 (range: 0.4-1.0)

### Network Validation

- **Modularity**: 0.4649 (excellent community structure)
- **Community purity**: 0.909 (clear gene-disease separation)
- **Network density**: 0.0203 (appropriate for biological networks)
- **Path efficiency**: Average path length of 3.05

### Statistical Significance

- **Evidence Index**: Mean 0.922 ± 0.225 (strong literature support)
- **Multi-gene diseases**: 38 diseases connected to multiple genes
- **Scale-free properties**: Partial power-law distribution (R² = 0.696)

## Key Insights

- Each gene forms its own community with a specific set of associated diseases
- Some diseases (notably Alzheimer's) are linked to multiple genes
- The graph structure reveals patterns in gene-disease relationships
- High-quality associations show strong literature evidence (30.3% high-confidence)
- Network exhibits efficient small-world properties with clear modularity

## Example Results

When analyzing five genes (APP, MAPT, PSEN1, SNCA, BACE1) associated with neurodegenerative diseases:

- Alzheimer's Disease is associated with all five genes
- Multiple disease connections exist between specific gene pairs
- Graph community detection identifies gene-specific clusters
- Network validation confirms biological significance with 0.4649 modularity
- 38 diseases show multi-gene associations indicating pathway convergence

For detailed documentation, code explanations, and analysis results, see [DOCUMENTATION.md](./DOCUMENTATION.md).

## Generated Outputs

The analysis generates several output files in the `output/` directory:

### Visualizations
- `gene_disease_graph.png`: Single gene (APP) network visualization
- `multi_gene_disease_graph.png`: Multi-gene network with 5 genes
- `gene_disease_communities.png`: Community structure visualization  
- `validation_metrics_summary.png`: Comprehensive validation metrics

### Data Files
- `gene_disease_graph.graphml`: Single gene network in GraphML format
- `multi_gene_disease_graph.graphml`: Multi-gene network in GraphML format
- `validation_results.json`: Detailed validation and precision metrics
- `summary_report.md`: Comprehensive analysis summary report

## Project Structure

```
Gene-Disease-Graph/
├── collectData.py              # Data collection from DisGeNET API
├── create_graph.py             # Single gene graph creation
├── multi_gene_graph.py         # Multi-gene graph creation  
├── analyze_communities.py      # Community detection analysis
├── advanced_validation.py      # Precision metrics and validation
├── generate_summary.py         # Comprehensive summary generator
├── README.md                   # Project overview (this file)
├── DOCUMENTATION.md            # Detailed documentation
├── .env                        # API key configuration
├── data/                       # Raw JSON data from DisGeNET
└── output/                     # Generated graphs, metrics, and reports
```

## Research Applications

This framework supports various research applications:

- **Drug Discovery**: Identify potential therapeutic targets and repurposing opportunities
- **Precision Medicine**: Develop gene-specific treatment strategies  
- **Biomarker Discovery**: Find common pathways for diagnostic markers
- **Pathway Analysis**: Understand genetic convergence in disease mechanisms
- **Clinical Translation**: Support evidence-based precision medicine approaches
