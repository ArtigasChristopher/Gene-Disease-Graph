import json
import networkx as nx
import matplotlib.pyplot as plt
import os
from collections import defaultdict

def load_data(file_path):
    """
    Load JSON data from a file
    """
    with open(file_path, 'r', encoding='utf-8') as f:
        return json.load(f)

def create_gene_disease_graph(data):
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

    # Iterate through gene-disease associations
    for assoc in data['payload']:
        # Check if record contains necessary information
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
    print(f"Complete records: {complete_records} out of {len(data['payload'])}")
    print(f"Total associations: {total_associations}")

    return G, gene_info, disease_info

def visualize_graph(G, title="Gene-Disease Graph", layout=None, labels=None, node_color=None):
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

    # Set graph layout
    if layout is None:
        layout = nx.spring_layout(G)

    # Set node labels
    if labels is None:
        labels = {}
        for node in G.nodes():
            if G.nodes[node].get('type') == 'gene':
                labels[node] = G.nodes[node].get('symbol', '')
            else:
                labels[node] = G.nodes[node].get('name', '')

    # Draw the graph
    nx.draw_networkx(G, pos=layout, node_color=node_color,
                     node_size=500, font_size=8, labels=labels,
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

    return metrics, node_types, top_degree_nodes

def main():
    # Ensure the output directory exists
    output_dir = os.path.join(os.path.dirname(__file__), 'output')
    os.makedirs(output_dir, exist_ok=True)

    # Load data
    data_path = os.path.join(os.path.dirname(__file__), 'data', 'summary_response.json')
    data = load_data(data_path)

    # Create the graph
    G, gene_info, disease_info = create_gene_disease_graph(data)

    # Analyze the graph
    analyze_graph(G)

    # Visualize the graph
    layout = visualize_graph(G)

    # Save the visualization
    plt.savefig(os.path.join(output_dir, 'gene_disease_graph.png'), dpi=300, bbox_inches='tight')
    print(f"Visualization saved to output/gene_disease_graph.png")

    # Save the graph
    save_graph(G, os.path.join(output_dir, 'gene_disease_graph.graphml'))

    # Show the graph
    plt.show()

if __name__ == "__main__":
    main()
