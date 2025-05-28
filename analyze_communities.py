import networkx as nx
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import os
import json
from collections import defaultdict
import community as community_louvain  # python-louvain package

def load_data(file_path):
    """
    Load JSON data from a file
    """
    with open(file_path, 'r', encoding='utf-8') as f:
        return json.load(f)

def load_graphml(file_path):
    """
    Load a graph from a GraphML file
    """
    G = nx.read_graphml(file_path)

    # Convert certain attributes to appropriate types
    for u, v, data in G.edges(data=True):
        if 'score' in data:
            try:
                data['score'] = float(data['score'])
            except (ValueError, TypeError):
                pass
        if 'ei' in data:
            try:
                data['ei'] = float(data['ei'])
            except (ValueError, TypeError):
                pass
        if 'num_pmids' in data:
            try:
                data['num_pmids'] = int(data['num_pmids'])
            except (ValueError, TypeError):
                pass

    return G

def detect_communities(G):
    """
    Detect communities in the graph using the Louvain algorithm
    """
    # Apply the Louvain community detection algorithm
    partition = community_louvain.best_partition(G)

    # Count the number of communities and calculate modularity
    num_communities = len(set(partition.values()))
    modularity = community_louvain.modularity(partition, G)

    print(f"Number of detected communities: {num_communities}")
    print(f"Modularity: {modularity:.4f}")

    return partition

def visualize_communities(G, partition, output_path=None):
    """
    Visualize the graph with different colors for each community
    """
    plt.figure(figsize=(15, 12))

    # Define a color for each community
    cmap = cm.get_cmap('tab20', max(partition.values()) + 1)

    # Define the position of the nodes using the spring_layout algorithm
    pos = nx.spring_layout(G, k=0.15, iterations=50, seed=42)

    # Draw the nodes
    for comm in set(partition.values()):
        list_nodes = [node for node in G.nodes() if partition[node] == comm]
        nx.draw_networkx_nodes(G, pos,
                              nodelist=list_nodes,
                              node_color=[cmap(comm)] * len(list_nodes),
                              node_size=[300 if G.nodes[node].get('type') == 'gene' else
                                        100 + G.degree[node] * 20 for node in list_nodes],
                              alpha=0.8)

    # Draw the edges
    nx.draw_networkx_edges(G, pos, width=0.5, alpha=0.5, edge_color='gray')

    # Add labels only for important nodes
    labels = {}
    for node in G.nodes():
        if G.nodes[node].get('type') == 'gene':
            labels[node] = G.nodes[node].get('symbol', '')
        elif G.degree[node] > 2:  # Show only diseases connected to multiple genes
            labels[node] = G.nodes[node].get('name', '')

    nx.draw_networkx_labels(G, pos, labels=labels, font_size=8, font_weight='bold')

    plt.title("Communities in the gene-disease graph", fontsize=16)
    plt.axis('off')

    if output_path:
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        print(f"Community visualization saved to {output_path}")

    return pos

def analyze_communities(G, partition):
    """
    Analyze the detected communities
    """
    # Group nodes by community
    communities = defaultdict(list)
    for node, comm_id in partition.items():
        communities[comm_id].append(node)

    # Analyze each community
    print("\nCommunity analysis:")

    for comm_id, nodes in sorted(communities.items(), key=lambda x: len(x[1]), reverse=True):
        genes = [node for node in nodes if G.nodes[node].get('type') == 'gene']
        diseases = [node for node in nodes if G.nodes[node].get('type') == 'disease']

        # Get the names/symbols of the genes and diseases
        gene_names = [G.nodes[g].get('symbol', '') for g in genes]
        disease_names = []
        for d in diseases:
            if len(disease_names) < 5:  # Limit to 5 disease names for readability
                disease_names.append(G.nodes[d].get('name', ''))

        print(f"\nCommunity {comm_id} ({len(nodes)} nodes: {len(genes)} genes, {len(diseases)} diseases)")
        print(f"  Genes: {', '.join(gene_names)}")
        print(f"  Main diseases: {', '.join(disease_names[:5])}{' ...' if len(disease_names) > 5 else ''}")

        # Identify diseases connected to multiple genes in this community
        if len(genes) > 1:
            shared_diseases = []
            for d in diseases:
                gene_connections = [g for g in G.neighbors(d) if g in genes]
                if len(gene_connections) > 1:
                    shared_diseases.append((d, gene_connections))

            if shared_diseases:
                print("  Diseases shared by multiple genes in this community:")
                for d, connected_genes in sorted(shared_diseases, key=lambda x: len(x[1]), reverse=True)[:3]:
                    connected_gene_names = [G.nodes[g].get('symbol', '') for g in connected_genes]
                    print(f"    {G.nodes[d].get('name', '')}: {', '.join(connected_gene_names)}")

def main():
    # Load the multi-gene graph
    graph_path = os.path.join(os.path.dirname(__file__), 'output', 'multi_gene_disease_graph.graphml')
    output_dir = os.path.join(os.path.dirname(__file__), 'output')

    if not os.path.exists(graph_path):
        print(f"The file {graph_path} does not exist. Please run multi_gene_graph.py first.")
        return

    try:
        # Load the graph
        print(f"Loading graph from {graph_path}")
        G = load_graphml(graph_path)
        print(f"Graph loaded with {G.number_of_nodes()} nodes and {G.number_of_edges()} edges")

        # Detect communities
        try:
            partition = detect_communities(G)

            # Visualize communities
            output_path = os.path.join(output_dir, 'gene_disease_communities.png')
            visualize_communities(G, partition, output_path)

            # Analyze communities
            analyze_communities(G, partition)

            # Show the graph
            plt.show()

        except ImportError:
            print("The python-louvain package is not installed. Installing...")
            print("Use the command: pip install python-louvain")

    except Exception as e:
        print(f"Error while analyzing the graph: {e}")

if __name__ == "__main__":
    main()
