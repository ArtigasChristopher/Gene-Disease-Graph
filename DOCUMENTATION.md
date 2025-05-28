# Gene-Disease Association Graph

## Project Overview

This project aims to retrieve gene-disease association data from the DisGeNET database via its API and structure this data as a network graph. The graph represents genes and diseases as nodes, with edges representing the associations between them. This visualization and analysis approach helps in understanding the complex relationships between genes and diseases, particularly focusing on neurodegenerative disorders.

## Table of Contents

1. [Prerequisites](#prerequisites)
2. [Project Structure](#project-structure)
3. [Data Collection](#data-collection)
4. [Graph Creation and Visualization](#graph-creation-and-visualization)
5. [Multi-Gene Analysis](#multi-gene-analysis)
6. [Community Detection](#community-detection)
7. [Results and Insights](#results-and-insights)
8. [Future Enhancements](#future-enhancements)

## Prerequisites

To run this project, you need:

1. Python 3.x
2. DisGeNET API key (stored in a `.env` file)
3. Required Python packages:
   - requests
   - networkx
   - matplotlib
   - python-dotenv
   - python-louvain (for community detection)
   - pandas (for data manipulation)

You can install the required packages using pip:

```bash
pip install requests networkx matplotlib python-dotenv python-louvain pandas
```

## Project Structure

The project is organized as follows:

```
Gene-Disease-Graph/
├── collectData.py            # Script to collect data from DisGeNET API
├── create_graph.py           # Script to create and visualize a graph for a single gene
├── multi_gene_graph.py       # Script to create a graph for multiple genes
├── analyze_communities.py    # Script to detect and analyze communities in the graph
├── README.md                 # Project overview
├── DOCUMENTATION.md          # Detailed project documentation (this file)
├── .env                      # Contains DisGeNET API key (not tracked in git)
├── data/                     # Directory containing raw data from DisGeNET
│   ├── summary_response.json # Initial data for APP gene
│   ├── gene_351_summary.json # Data for APP gene
│   ├── gene_4137_summary.json # Data for MAPT gene
│   ├── gene_5663_summary.json # Data for PSEN1 gene
│   ├── gene_6622_summary.json # Data for SNCA gene
│   └── gene_23621_summary.json # Data for BACE1 gene
└── output/                   # Directory containing generated outputs
    ├── gene_disease_graph.png # Visualization of single gene graph
    ├── gene_disease_graph.graphml # GraphML file for single gene
    ├── multi_gene_disease_graph.png # Visualization of multi-gene graph
    ├── multi_gene_disease_graph.graphml # GraphML file for multi-gene graph
    └── gene_disease_communities.png # Visualization of communities
```

## Data Collection

The first step in the project is to collect gene-disease association data from the DisGeNET database using their API.

### API Key Setup

The DisGeNET API key is stored in a `.env` file at the project root:

```
DISGENET_API_KEY = "your_api_key_here"
```

### Data Collection Script

The `collectData.py` script fetches gene-disease association data for a specific gene (initially the APP gene, NCBI ID: 351) and saves it as JSON:

```python
import requests
import json
import time
from dotenv import load_dotenv
import os

def load_api_key():
    load_dotenv()
    api_key = os.getenv("DISGENET_API_KEY")
    if api_key is None:
        raise ValueError("API key not found. Please set the DISGENET_API_KEY environment variable in your .env file.")
    return api_key

def get_gda_summary(api_key, params, max_retries=1):
    headers = {
        "Authorization": api_key,
        "Accept": "application/json"
    }
    url = "https://api.disgenet.com/api/v1/gda/summary"
    for attempt in range(max_retries + 1):
        try:
            response = requests.get(
                url,
                params=params,
                headers=headers,
                timeout=10
            )
        except requests.RequestException as e:
            print(f"Request failed: {e}")
            if attempt < max_retries:
                time.sleep(2)
                continue
            else:
                raise

        if response.status_code == 429 and attempt < max_retries:
            retry_after = int(response.headers.get("x-rate-limit-retry-after-seconds", 60))
            print(f"Rate limit reached, waiting {retry_after}s...")
            time.sleep(retry_after)
            continue
        response.raise_for_status()
        return response.json()
    raise Exception("Failed to get a valid response after retries.")

def save_json(data, filename):
    with open(filename, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=2, ensure_ascii=False)
    print(f"Summary response saved to {filename}")

def main():
    api_key = load_api_key()
    params = {
        "gene_ncbi_id": "351",  # APP gene
        "page_number": "0"
    }
    response_parsed = get_gda_summary(api_key, params, max_retries=1)
    save_json(response_parsed, "data/summary_response.json")

if __name__ == "__main__":
    main()
```

## Graph Creation and Visualization

After collecting the data, we create a graph representation of the gene-disease associations using NetworkX.

### Single Gene Graph

The `create_graph.py` script processes the collected data for a single gene and creates a graph:

```python
import json
import networkx as nx
import matplotlib.pyplot as plt
import os
from collections import defaultdict

# Function definitions for loading data, creating graph, etc.
# ...

def main():
    # Load data
    data_path = os.path.join(os.path.dirname(__file__), 'data', 'summary_response.json')
    data = load_data(data_path)

    # Create graph
    G, gene_info, disease_info = create_gene_disease_graph(data)

    # Analyze graph
    analyze_graph(G)

    # Visualize graph
    layout = visualize_graph(G)

    # Save visualization and graph
    plt.savefig(os.path.join(output_dir, 'gene_disease_graph.png'), dpi=300, bbox_inches='tight')
    save_graph(G, os.path.join(output_dir, 'gene_disease_graph.graphml'))

    # Display graph
    plt.show()
```

### Example Output for Single Gene

When running the single gene graph script for the APP gene (ID: 351), we get:

```
Graph created with 58 nodes and 57 edges
Complete records: 57 out of 57
Total associations: 57
Graph metrics:
Number of nodes: 58
Number of edges: 57
Graph density: 0.034483
Average degree: 1.97

Node types:
gene: 1
disease: 57

Most connected nodes:
Gene APP: 57 connections
Disease Alzheimer's Disease: 1 connections
Disease Impaired cognition: 1 connections
Disease CEREBRAL AMYLOID ANGIOPATHY, APP-RELATED: 1 connections
...
```

The visualization shows a star-like structure with the APP gene in the center connected to various diseases.

## Multi-Gene Analysis

To gain deeper insights, we extend the analysis to multiple genes related to neurodegenerative diseases.

### Multi-Gene Graph Script

The `multi_gene_graph.py` script fetches data for multiple genes and creates a more complex graph:

```python
import requests
import json
import time
import os
import networkx as nx
import matplotlib.pyplot as plt
from dotenv import load_dotenv
from collections import defaultdict
import random

# Function definitions...

def main():
    # List of gene IDs to fetch
    # APP (351) + some genes associated with neurodegenerative diseases
    gene_ids = [351, 4137, 5663, 6622, 23621]  # APP, MAPT, PSEN1, SNCA, BACE1

    # Fetch data for all genes
    data_list = fetch_multiple_genes(gene_ids, data_dir)

    # Create graph
    G, gene_info, disease_info = create_gene_disease_graph(data_list)

    # Analyze, visualize, and save graph
    # ...
```

### Example Output for Multi-Gene Analysis

The multi-gene analysis provides richer insights:

```
Graph created with 135 nodes and 184 edges
Complete records: 184
Total associations: 184
Graph metrics:
Number of nodes: 135
Number of edges: 184
Graph density: 0.020343
Average degree: 2.73

Node types:
gene: 5
disease: 130

Most connected nodes:
Gene APP: 57 connections
Gene PSEN1: 48 connections
Gene MAPT: 46 connections
Gene SNCA: 31 connections
Disease Alzheimer's Disease: 5 connections
...

Diseases shared by multiple genes:
Alzheimer's Disease: APP, MAPT, PSEN1, SNCA, BACE1
Neurodegenerative Disorders: APP, MAPT, PSEN1, SNCA
Degenerative disorder (disorder): APP, MAPT, PSEN1, SNCA
Memory Disorders: APP, MAPT, PSEN1
Nerve Degeneration: APP, PSEN1, SNCA
```

## Community Detection

To identify clusters in the graph, we perform community detection analysis.

### Community Analysis Script

The `analyze_communities.py` script applies the Louvain community detection algorithm:

```python
import networkx as nx
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import os
import json
from collections import defaultdict
import community as community_louvain  # python-louvain package

# Function definitions...

def main():
    # Load multi-gene graph
    graph_path = os.path.join(os.path.dirname(__file__), 'output', 'multi_gene_disease_graph.graphml')

    # Load graph
    G = load_graphml(graph_path)

    # Detect communities
    partition = detect_communities(G)

    # Visualize communities
    output_path = os.path.join(output_dir, 'gene_disease_communities.png')
    visualize_communities(G, partition, output_path)

    # Analyze communities
    analyze_communities(G, partition)
```

### Example Output for Community Detection

The community detection analysis found:

```
Number of detected communities: 5
Modularity: 0.4649

Community analysis:

Community 0 (36 nodes: 1 genes, 35 diseases)
  Genes: APP
  Main diseases: Impaired cognition, CEREBRAL AMYLOID ANGIOPATHY, APP-RELATED, Fragile X Syndrome, Cerebral Amyloid Angiopathy, Hereditary, Cerebral hemorrhage with amyloidosis, hereditary, Dutch type

Community 3 (33 nodes: 1 genes, 32 diseases)
  Genes: PSEN1
  Main diseases: familial Alzheimer disease, Plaque, Amyloid, Amyloidosis, Depressive disorder, Mental Depression
...
```

This reveals that each gene forms its own community with its associated diseases.

## Results and Insights

From our analysis, we can draw several key insights:

1. **Gene-Specific Disease Profiles**: Each gene has a specific set of diseases it's associated with, forming distinct communities in the graph.

2. **Shared Diseases**: Some diseases, particularly Alzheimer's Disease, are associated with multiple genes (APP, MAPT, PSEN1, SNCA, BACE1), suggesting they involve multiple genetic pathways.

3. **Disease Connectivity**: While most diseases are associated with only one gene in our dataset, a few (like Alzheimer's Disease and Neurodegenerative Disorders) are connected to multiple genes, indicating their complex genetic basis.

4. **Community Structure**: The network naturally divides into communities that align with individual genes, suggesting that genetic factors create natural clustering in disease associations.

## Future Enhancements

Several further developments could enhance this project:

1. **Expanded Gene Set**: Include more genes associated with neurodegenerative and other diseases to build a more comprehensive network.

2. **Additional Data Sources**: Incorporate data from other sources such as gene expression databases, protein-protein interaction networks, or pathway databases.

3. **Temporal Analysis**: Analyze how gene-disease associations have evolved over time using the publication year data.

4. **Path Analysis**: Explore paths between genes through shared diseases to identify potential biological connections.

5. **Interactive Visualization**: Create interactive web-based visualizations to explore the graph more effectively.

6. **Machine Learning Integration**: Apply machine learning algorithms to predict new gene-disease associations based on graph patterns.

---

This documentation provides a comprehensive guide to understanding and using the Gene-Disease Association Graph project. The code, visualizations, and analyses presented here offer insights into the complex relationships between genes and diseases, with a particular focus on neurodegenerative disorders.
