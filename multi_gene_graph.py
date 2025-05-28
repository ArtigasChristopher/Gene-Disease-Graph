import requests
import json
import time
import os
import networkx as nx
import matplotlib.pyplot as plt
from dotenv import load_dotenv
from collections import defaultdict
import random

def load_api_key():
    """
    Load API key from .env file
    """
    load_dotenv()
    api_key = os.getenv("DISGENET_API_KEY")
    if api_key is None:
        raise ValueError("API key not found. Please set the DISGENET_API_KEY environment variable in your .env file.")
    return api_key

def get_gda_summary(api_key, params, max_retries=3):
    """
    Retrieve gene-disease association data via DisGeNET API
    """
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
    """
    Save data as JSON format
    """
    with open(filename, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=2, ensure_ascii=False)
    print(f"Data saved to {filename}")

def load_data(file_path):
    """
    Load JSON data from a file
    """
    with open(file_path, 'r', encoding='utf-8') as f:
        return json.load(f)

def create_gene_disease_graph(data_list):
    """
    Create a graph from gene-disease association data
    """
    # Initialize the graph
    G = nx.Graph()

    # Dictionaries to store node information
    gene_info = {}
    disease_info = {}

    # Counters for statistics
    total_associations = 0
    complete_records = 0

    # Iterate over gene-disease associations for each dataset
    for data in data_list:
        for assoc in data['payload']:
            # Check if the record contains the necessary information
            if 'symbolOfGene' not in assoc or 'diseaseName' not in assoc:
                continue

            complete_records += 1

            # Extract gene and disease information
            gene_id = assoc['geneNcbiID']
            gene_symbol = assoc['symbolOfGene']
            disease_id = assoc.get('diseaseUMLSCUI', '')
            disease_name = assoc['diseaseName']

            # Add nodes for genes and diseases if they don't already exist
            if gene_id not in gene_info:
                G.add_node(f"gene_{gene_id}",
                          type='gene',
                          symbol=gene_symbol,
                          id=gene_id)
                gene_info[gene_id] = {
                    'symbol': gene_symbol,
                    'dsi': assoc.get('geneDSI', None),
                    'dpi': assoc.get('geneDPI', None),
                    'pli': assoc.get('genepLI', None)
                }

            if disease_id and disease_id not in disease_info:
                disease_classes = []
                for class_type in ['diseaseClasses_MSH', 'diseaseClasses_DO', 'diseaseClasses_HPO', 'diseaseClasses_UMLS_ST']:
                    if class_type in assoc and assoc[class_type]:
                        disease_classes.extend(assoc[class_type])

                G.add_node(f"disease_{disease_id}",
                          type='disease',
                          name=disease_name,
                          id=disease_id,
                          classes=disease_classes)
                disease_info[disease_id] = {
                    'name': disease_name,
                    'type': assoc.get('diseaseType', ''),
                    'classes': disease_classes
                }

            # Add an edge for the gene-disease association
            if disease_id:  # Ensure a disease ID is available
                G.add_edge(f"gene_{gene_id}", f"disease_{disease_id}",
                          score=assoc.get('score', 0),
                          ei=assoc.get('ei', 0),
                          year_initial=assoc.get('yearInitial', None),
                          year_final=assoc.get('yearFinal', None),
                          num_pmids=assoc.get('numPMIDs', 0))
                total_associations += 1

    print(f"Graph created with {G.number_of_nodes()} nodes and {G.number_of_edges()} edges")
    print(f"Complete records: {complete_records}")
    print(f"Total associations: {total_associations}")

    return G, gene_info, disease_info

def visualize_graph(G, title="Gene-Disease Graph", layout=None, labels=None, node_color=None, node_size=None):
    """
    Visualize the graph with formatting options
    """
    plt.figure(figsize=(12, 10))

    # Set node colors based on their type
    if node_color is None:
        node_color = []
        for node in G.nodes():
            if G.nodes[node].get('type') == 'gene':
                node_color.append('skyblue')
            else:
                node_color.append('lightcoral')

    # Set node size based on their degree
    if node_size is None:
        node_size = []
        for node in G.nodes():
            degree = G.degree[node]
            if G.nodes[node].get('type') == 'gene':
                node_size.append(300 + degree * 5)  # Base size plus larger for genes
            else:
                node_size.append(100 + degree * 20)  # Size proportional to connections

    # Set graph layout
    if layout is None:
        layout = nx.spring_layout(G, k=0.15, iterations=50)

    # Set node labels
    if labels is None:
        labels = {}
        for node in G.nodes():
            if G.nodes[node].get('type') == 'gene':
                labels[node] = G.nodes[node].get('symbol', '')
            else:
                # To avoid label overcrowding, only show labels for important diseases
                if G.degree[node] > 1:
                    labels[node] = G.nodes[node].get('name', '')

    # Draw the graph
    nx.draw_networkx(G, pos=layout, node_color=node_color,
                     node_size=node_size, font_size=8, labels=labels,
                     edge_color='gray', width=0.5, alpha=0.7)

    plt.title(title)
    plt.axis('off')

    return layout

def save_graph(G, file_path):
    """
    Save the graph in GraphML format
    """
    # Create a copy of the graph for saving
    G_save = nx.Graph()

    # Copy nodes, converting lists to strings and handling None values
    for node, attrs in G.nodes(data=True):
        node_attrs = {}
        for key, value in attrs.items():
            if value is None:
                node_attrs[key] = "N/A"  # Replace None with a string
            elif isinstance(value, list):
                node_attrs[key] = ';'.join(str(v) for v in value)
            else:
                node_attrs[key] = value
        G_save.add_node(node, **node_attrs)

    # Copy edges
    for u, v, attrs in G.edges(data=True):
        edge_attrs = {}
        for key, value in attrs.items():
            if value is None:
                edge_attrs[key] = "N/A"  # Replace None with a string
            elif isinstance(value, list):
                edge_attrs[key] = ';'.join(str(v) for v in value)
            else:
                edge_attrs[key] = value
        G_save.add_edge(u, v, **edge_attrs)

    nx.write_graphml(G_save, file_path)
    print(f"Graph saved to {file_path}")

def analyze_graph(G):
    """
    Analyze the graph to extract metrics and statistics
    """
    # Calculate basic metrics
    metrics = {
        'nodes': G.number_of_nodes(),
        'edges': G.number_of_edges(),
        'density': nx.density(G),
        'avg_degree': sum(dict(G.degree()).values()) / G.number_of_nodes()
    }

    # Count node types
    node_types = defaultdict(int)
    for node in G.nodes():
        node_type = G.nodes[node].get('type', 'unknown')
        node_types[node_type] += 1

    # Identify the most connected nodes
    degrees = dict(G.degree())
    top_degree_nodes = sorted(degrees.items(), key=lambda x: x[1], reverse=True)[:10]

    print("Graph metrics:")
    print(f"Number of nodes: {metrics['nodes']}")
    print(f"Number of edges: {metrics['edges']}")
    print(f"Graph density: {metrics['density']:.6f}")
    print(f"Average degree: {metrics['avg_degree']:.2f}")

    print("\nNode types:")
    for node_type, count in node_types.items():
        print(f"{node_type}: {count}")

    print("\nMost connected nodes:")
    for node, degree in top_degree_nodes:
        node_info = G.nodes[node]
        if node_info.get('type') == 'gene':
            print(f"Gene {node_info.get('symbol', '')}: {degree} connections")
        else:
            print(f"Disease {node_info.get('name', '')}: {degree} connections")

    # Identify diseases common to multiple genes
    if node_types['gene'] > 1:
        disease_genes = defaultdict(list)
        for disease_node in [n for n in G.nodes() if G.nodes[n].get('type') == 'disease']:
            for neighbor in G.neighbors(disease_node):
                if G.nodes[neighbor].get('type') == 'gene':
                    disease_genes[disease_node].append(neighbor)

        shared_diseases = {d: genes for d, genes in disease_genes.items() if len(genes) > 1}

        if shared_diseases:
            print("\nDiseases shared by multiple genes:")
            for disease, genes in sorted(shared_diseases.items(), key=lambda x: len(x[1]), reverse=True)[:5]:
                gene_symbols = [G.nodes[g].get('symbol', '') for g in genes]
                print(f"{G.nodes[disease].get('name', '')}: {', '.join(gene_symbols)}")

    return metrics, node_types, top_degree_nodes

def fetch_multiple_genes(gene_ids, output_dir="data"):
    """
    Retrieve association data for multiple genes
    """
    api_key = load_api_key()
    data_list = []

    for gene_id in gene_ids:
        output_file = os.path.join(output_dir, f"gene_{gene_id}_summary.json")

        # Check if data is already available
        if os.path.exists(output_file):
            print(f"Loading existing data for gene {gene_id}")
            data = load_data(output_file)
        else:
            print(f"Fetching data for gene {gene_id}")
            params = {
                "gene_ncbi_id": str(gene_id),
                "page_number": "0"
            }
            try:
                data = get_gda_summary(api_key, params, max_retries=3)
                save_json(data, output_file)
            except Exception as e:
                print(f"Error fetching data for gene {gene_id}: {e}")
                continue

        data_list.append(data)
        # Add a delay to respect API limits
        if len(gene_ids) > 1:
            time.sleep(1)

    return data_list

def main():
    # List of gene IDs to retrieve
    # APP (351) + some genes associated with neurodegenerative diseases
    gene_ids = [351, 4137, 5663, 6622, 23621]  # APP, MAPT, PSEN1, SNCA, BACE1

    # Ensure directories exist
    data_dir = os.path.join(os.path.dirname(__file__), 'data')
    output_dir = os.path.join(os.path.dirname(__file__), 'output')
    os.makedirs(data_dir, exist_ok=True)
    os.makedirs(output_dir, exist_ok=True)

    # Retrieve data for all genes
    data_list = fetch_multiple_genes(gene_ids, data_dir)

    if not data_list:
        print("No data retrieved.")
        return

    # Create the graph
    G, gene_info, disease_info = create_gene_disease_graph(data_list)

    # If the graph is empty, stop
    if G.number_of_nodes() == 0:
        print("The graph is empty.")
        return

    # Analyze the graph
    analyze_graph(G)

    # Visualize the graph
    layout = visualize_graph(G, title=f"Graph of gene-disease associations ({len(gene_info)} genes)")

    # Save the visualization
    plt.savefig(os.path.join(output_dir, 'multi_gene_disease_graph.png'), dpi=300, bbox_inches='tight')
    print(f"Visualization saved to output/multi_gene_disease_graph.png")

    # Save the graph
    save_graph(G, os.path.join(output_dir, 'multi_gene_disease_graph.graphml'))

    # Show the graph
    plt.show()

if __name__ == "__main__":
    main()
