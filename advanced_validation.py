#!/usr/bin/env python3
"""
Advanced Validation and Precision Metrics for Gene-Disease Network Analysis

This module provides comprehensive validation metrics and statistical analysis
for gene-disease association networks, including precision, recall, F1-score,
and network topology validation.
"""

import networkx as nx
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from sklearn.metrics import precision_score, recall_score, f1_score, roc_auc_score, confusion_matrix
from sklearn.model_selection import cross_val_score
from sklearn.preprocessing import LabelEncoder
import json
import os
from collections import defaultdict, Counter
import warnings
warnings.filterwarnings('ignore')

class GeneDiseasePrecisionValidator:
    """
    A comprehensive validator for gene-disease association networks
    """

    def __init__(self, graph_path, data_dir):
        """
        Initialize the validator with the graph and data
        """
        self.graph_path = graph_path
        self.data_dir = data_dir
        self.graph = self.load_graph()
        self.raw_data = self.load_raw_data()
        self.associations_df = self.create_associations_dataframe()
        self.validation_results = {}

    def load_graph(self):
        """Load the multi-gene graph"""
        G = nx.read_graphml(self.graph_path)

        # Convert attributes to appropriate types
        for u, v, data in G.edges(data=True):
            if 'score' in data:
                try:
                    data['score'] = float(data['score'])
                except (ValueError, TypeError):
                    data['score'] = 0.0
            if 'ei' in data:
                try:
                    data['ei'] = float(data['ei'])
                except (ValueError, TypeError):
                    data['ei'] = 0.0
            if 'num_pmids' in data:
                try:
                    data['num_pmids'] = int(data['num_pmids'])
                except (ValueError, TypeError):
                    data['num_pmids'] = 0

        print(f"Loaded graph with {G.number_of_nodes()} nodes and {G.number_of_edges()} edges")
        return G

    def load_raw_data(self):
        """Load all raw JSON data files"""
        raw_data = []
        json_files = [f for f in os.listdir(self.data_dir) if f.endswith('.json')]

        for file in json_files:
            file_path = os.path.join(self.data_dir, file)
            with open(file_path, 'r', encoding='utf-8') as f:
                data = json.load(f)
                raw_data.extend(data.get('payload', []))

        print(f"Loaded {len(raw_data)} raw associations from {len(json_files)} files")
        return raw_data

    def create_associations_dataframe(self):
        """Create a comprehensive DataFrame of all associations"""
        associations = []

        for assoc in self.raw_data:
            if 'symbolOfGene' in assoc and 'diseaseName' in assoc:
                associations.append({
                    'gene_id': assoc.get('geneNcbiID'),
                    'gene_symbol': assoc.get('symbolOfGene'),
                    'disease_id': assoc.get('diseaseUMLSCUI', ''),
                    'disease_name': assoc.get('diseaseName'),
                    'score': float(assoc.get('score', 0) or 0),
                    'ei': float(assoc.get('ei', 0) or 0),
                    'num_pmids': int(assoc.get('numPMIDs', 0) or 0),
                    'year_initial': assoc.get('yearInitial'),
                    'year_final': assoc.get('yearFinal'),
                    'disease_type': assoc.get('diseaseType', ''),
                    'gene_dsi': assoc.get('geneDSI'),
                    'gene_dpi': assoc.get('geneDPI'),
                    'gene_pli': assoc.get('genepLI')
                })

        df = pd.DataFrame(associations)
        print(f"Created associations DataFrame with {len(df)} records")
        return df

    def calculate_association_quality_metrics(self):
        """Calculate quality metrics for gene-disease associations"""
        print("\n=== Association Quality Metrics ===")

        metrics = {}

        # Score distribution analysis
        score_stats = self.associations_df['score'].describe()
        metrics['score_statistics'] = score_stats.to_dict()

        # Evidence Index (EI) analysis
        ei_stats = self.associations_df['ei'].describe()
        metrics['ei_statistics'] = ei_stats.to_dict()

        # PMID count analysis
        pmid_stats = self.associations_df['num_pmids'].describe()
        metrics['pmid_statistics'] = pmid_stats.to_dict()

        print(f"Score range: {score_stats['min']:.3f} - {score_stats['max']:.3f}")
        print(f"Mean score: {score_stats['mean']:.3f} ± {score_stats['std']:.3f}")
        print(f"Evidence Index range: {ei_stats['min']:.3f} - {ei_stats['max']:.3f}")
        print(f"Mean PMID count: {pmid_stats['mean']:.1f}")

        # High-confidence associations (top quartile)
        high_confidence_threshold = score_stats['75%']
        high_confidence_count = (self.associations_df['score'] >= high_confidence_threshold).sum()
        metrics['high_confidence_associations'] = high_confidence_count
        metrics['high_confidence_percentage'] = (high_confidence_count / len(self.associations_df)) * 100

        print(f"High-confidence associations (score ≥ {high_confidence_threshold:.3f}): {high_confidence_count} ({metrics['high_confidence_percentage']:.1f}%)")

        return metrics

    def validate_network_topology(self):
        """Validate the network topology and structure"""
        print("\n=== Network Topology Validation ===")

        topology_metrics = {}

        # Basic network metrics
        topology_metrics['nodes'] = self.graph.number_of_nodes()
        topology_metrics['edges'] = self.graph.number_of_edges()
        topology_metrics['density'] = nx.density(self.graph)
        topology_metrics['avg_clustering'] = nx.average_clustering(self.graph)

        # Degree distribution analysis
        degrees = [d for n, d in self.graph.degree()]
        topology_metrics['degree_mean'] = np.mean(degrees)
        topology_metrics['degree_std'] = np.std(degrees)
        topology_metrics['degree_max'] = max(degrees)

        # Component analysis
        if nx.is_connected(self.graph):
            topology_metrics['connected_components'] = 1
            topology_metrics['diameter'] = nx.diameter(self.graph)
            topology_metrics['avg_path_length'] = nx.average_shortest_path_length(self.graph)
        else:
            components = list(nx.connected_components(self.graph))
            topology_metrics['connected_components'] = len(components)
            largest_component = max(components, key=len)
            subgraph = self.graph.subgraph(largest_component)
            topology_metrics['largest_component_size'] = len(largest_component)
            topology_metrics['diameter'] = nx.diameter(subgraph)
            topology_metrics['avg_path_length'] = nx.average_shortest_path_length(subgraph)

        # Scale-free network test (power-law degree distribution)
        degree_counts = Counter(degrees)
        x = list(degree_counts.keys())
        y = list(degree_counts.values())

        if len(x) > 3:  # Need enough points for regression
            # Log-log regression to test power-law
            log_x = np.log(x)
            log_y = np.log(y)
            slope, intercept, r_value, p_value, std_err = stats.linregress(log_x, log_y)

            topology_metrics['power_law_exponent'] = -slope
            topology_metrics['power_law_r_squared'] = r_value**2
            topology_metrics['power_law_p_value'] = p_value

            is_scale_free = (r_value**2 > 0.8 and p_value < 0.05 and -3 < slope < -1)
            topology_metrics['is_scale_free'] = is_scale_free

        print(f"Network density: {topology_metrics['density']:.4f}")
        print(f"Average clustering coefficient: {topology_metrics['avg_clustering']:.4f}")
        print(f"Average degree: {topology_metrics['degree_mean']:.2f} ± {topology_metrics['degree_std']:.2f}")
        print(f"Connected components: {topology_metrics['connected_components']}")
        print(f"Network diameter: {topology_metrics['diameter']}")
        print(f"Average path length: {topology_metrics['avg_path_length']:.2f}")

        if 'is_scale_free' in topology_metrics:
            print(f"Scale-free network: {'Yes' if topology_metrics['is_scale_free'] else 'No'} (R² = {topology_metrics['power_law_r_squared']:.3f})")

        return topology_metrics

    def analyze_gene_disease_specificity(self):
        """Analyze specificity of gene-disease associations"""
        print("\n=== Gene-Disease Specificity Analysis ===")

        specificity_metrics = {}

        # Gene specificity (how many diseases per gene)
        gene_disease_counts = self.associations_df.groupby('gene_symbol')['disease_name'].nunique()
        specificity_metrics['diseases_per_gene'] = gene_disease_counts.to_dict()
        specificity_metrics['avg_diseases_per_gene'] = gene_disease_counts.mean()
        specificity_metrics['std_diseases_per_gene'] = gene_disease_counts.std()

        # Disease specificity (how many genes per disease)
        disease_gene_counts = self.associations_df.groupby('disease_name')['gene_symbol'].nunique()
        specificity_metrics['genes_per_disease'] = disease_gene_counts.to_dict()
        specificity_metrics['avg_genes_per_disease'] = disease_gene_counts.mean()
        specificity_metrics['std_genes_per_disease'] = disease_gene_counts.std()

        # Multi-gene diseases (potential hubs)
        multi_gene_diseases = disease_gene_counts[disease_gene_counts > 1]
        specificity_metrics['multi_gene_diseases_count'] = len(multi_gene_diseases)
        specificity_metrics['multi_gene_diseases'] = multi_gene_diseases.to_dict()

        print(f"Average diseases per gene: {specificity_metrics['avg_diseases_per_gene']:.1f} ± {specificity_metrics['std_diseases_per_gene']:.1f}")
        print(f"Average genes per disease: {specificity_metrics['avg_genes_per_disease']:.1f} ± {specificity_metrics['std_genes_per_disease']:.1f}")
        print(f"Diseases associated with multiple genes: {specificity_metrics['multi_gene_diseases_count']}")

        if specificity_metrics['multi_gene_diseases_count'] > 0:
            print("Top multi-gene diseases:")
            for disease, count in sorted(multi_gene_diseases.items(), key=lambda x: x[1], reverse=True)[:5]:
                print(f"  {disease}: {count} genes")

        return specificity_metrics

    def calculate_precision_recall_metrics(self):
        """Calculate precision and recall based on score thresholds"""
        print("\n=== Precision-Recall Analysis ===")

        # Define different score thresholds
        thresholds = np.percentile(self.associations_df['score'], [50, 60, 70, 80, 90, 95])

        precision_recall_metrics = {}

        for threshold in thresholds:
            # Binary classification: high score (1) vs low score (0)
            y_true = (self.associations_df['score'] >= threshold).astype(int)
            y_pred = (self.associations_df['ei'] >= self.associations_df['ei'].median()).astype(int)

            if len(np.unique(y_true)) > 1:  # Check if we have both classes
                precision = precision_score(y_true, y_pred, zero_division=0)
                recall = recall_score(y_true, y_pred, zero_division=0)
                f1 = f1_score(y_true, y_pred, zero_division=0)

                precision_recall_metrics[f'threshold_{threshold:.2f}'] = {
                    'precision': precision,
                    'recall': recall,
                    'f1_score': f1,
                    'support': y_true.sum()
                }

        # Overall quality assessment based on multiple evidence sources
        # High-quality associations: high score AND multiple PMIDs
        high_quality = ((self.associations_df['score'] >= self.associations_df['score'].quantile(0.75)) &
                       (self.associations_df['num_pmids'] >= 2))

        medium_quality = ((self.associations_df['score'] >= self.associations_df['score'].median()) &
                         (self.associations_df['num_pmids'] >= 1))

        precision_recall_metrics['quality_distribution'] = {
            'high_quality': high_quality.sum(),
            'medium_quality': medium_quality.sum(),
            'total_associations': len(self.associations_df),
            'high_quality_percentage': (high_quality.sum() / len(self.associations_df)) * 100,
            'medium_quality_percentage': (medium_quality.sum() / len(self.associations_df)) * 100
        }

        print(f"High-quality associations: {high_quality.sum()} ({precision_recall_metrics['quality_distribution']['high_quality_percentage']:.1f}%)")
        print(f"Medium-quality associations: {medium_quality.sum()} ({precision_recall_metrics['quality_distribution']['medium_quality_percentage']:.1f}%)")

        return precision_recall_metrics

    def validate_community_structure(self):
        """Validate the community detection results"""
        print("\n=== Community Structure Validation ===")

        try:
            import community as community_louvain
        except ImportError:
            print("Community detection package not available")
            return {}

        # Detect communities
        partition = community_louvain.best_partition(self.graph)
        modularity = community_louvain.modularity(partition, self.graph)

        community_metrics = {
            'modularity': modularity,
            'num_communities': len(set(partition.values()))
        }

        # Analyze community composition
        communities = defaultdict(list)
        for node, comm_id in partition.items():
            communities[comm_id].append(node)

        # Calculate community purity (gene vs disease separation)
        community_purities = []
        for comm_id, nodes in communities.items():
            gene_count = sum(1 for node in nodes if self.graph.nodes[node].get('type') == 'gene')
            disease_count = len(nodes) - gene_count

            if len(nodes) > 0:
                # Purity = fraction of the majority type
                purity = max(gene_count, disease_count) / len(nodes)
                community_purities.append(purity)

        community_metrics['avg_community_purity'] = np.mean(community_purities)
        community_metrics['community_purities'] = community_purities

        # Silhouette-like analysis for communities
        # Calculate intra-community vs inter-community edge density
        intra_densities = []
        for comm_id, nodes in communities.items():
            subgraph = self.graph.subgraph(nodes)
            if len(nodes) > 1:
                intra_densities.append(nx.density(subgraph))

        community_metrics['avg_intra_community_density'] = np.mean(intra_densities) if intra_densities else 0

        print(f"Modularity: {modularity:.4f}")
        print(f"Number of communities: {community_metrics['num_communities']}")
        print(f"Average community purity: {community_metrics['avg_community_purity']:.3f}")
        print(f"Average intra-community density: {community_metrics['avg_intra_community_density']:.4f}")

        return community_metrics

    def generate_comprehensive_report(self, output_dir):
        """Generate a comprehensive validation report"""
        print("\n=== Generating Comprehensive Validation Report ===")

        # Run all validation analyses
        quality_metrics = self.calculate_association_quality_metrics()
        topology_metrics = self.validate_network_topology()
        specificity_metrics = self.analyze_gene_disease_specificity()
        precision_recall_metrics = self.calculate_precision_recall_metrics()
        community_metrics = self.validate_community_structure()

        # Compile all results
        self.validation_results = {
            'quality_metrics': quality_metrics,
            'topology_metrics': topology_metrics,
            'specificity_metrics': specificity_metrics,
            'precision_recall_metrics': precision_recall_metrics,
            'community_metrics': community_metrics
        }

        # Save detailed results
        results_file = os.path.join(output_dir, 'validation_results.json')
        with open(results_file, 'w') as f:
            json.dump(self.validation_results, f, indent=2, default=str)

        print(f"Detailed validation results saved to {results_file}")

        # Create visualizations
        self.create_validation_visualizations(output_dir)

        return self.validation_results

    def create_validation_visualizations(self, output_dir):
        """Create comprehensive validation visualizations"""
        # Set style
        plt.style.use('seaborn-v0_8')

        # Create figure with subplots
        fig, axes = plt.subplots(2, 3, figsize=(18, 12))
        fig.suptitle('Gene-Disease Network Validation Metrics', fontsize=16, fontweight='bold')

        # 1. Score distribution
        axes[0, 0].hist(self.associations_df['score'], bins=30, alpha=0.7, color='skyblue', edgecolor='black')
        axes[0, 0].set_xlabel('Association Score')
        axes[0, 0].set_ylabel('Frequency')
        axes[0, 0].set_title('Distribution of Association Scores')
        axes[0, 0].axvline(self.associations_df['score'].median(), color='red', linestyle='--', label='Median')
        axes[0, 0].axvline(self.associations_df['score'].quantile(0.75), color='orange', linestyle='--', label='Q3')
        axes[0, 0].legend()

        # 2. PMID count distribution
        pmid_counts = self.associations_df['num_pmids'].value_counts().sort_index()
        axes[0, 1].bar(pmid_counts.index[:10], pmid_counts.values[:10], alpha=0.7, color='lightcoral')
        axes[0, 1].set_xlabel('Number of PMIDs')
        axes[0, 1].set_ylabel('Frequency')
        axes[0, 1].set_title('Distribution of Evidence (PMID Count)')

        # 3. Gene-Disease specificity
        gene_disease_counts = self.associations_df.groupby('gene_symbol')['disease_name'].nunique()
        axes[0, 2].bar(range(len(gene_disease_counts)), sorted(gene_disease_counts.values, reverse=True),
                      alpha=0.7, color='lightgreen')
        axes[0, 2].set_xlabel('Genes (ranked)')
        axes[0, 2].set_ylabel('Number of Associated Diseases')
        axes[0, 2].set_title('Gene-Disease Association Specificity')

        # 4. Score vs Evidence correlation
        axes[1, 0].scatter(self.associations_df['score'], self.associations_df['ei'],
                          alpha=0.6, color='purple', s=20)
        axes[1, 0].set_xlabel('Association Score')
        axes[1, 0].set_ylabel('Evidence Index')
        axes[1, 0].set_title('Score vs Evidence Index Correlation')

        # Calculate correlation
        correlation = self.associations_df['score'].corr(self.associations_df['ei'])
        axes[1, 0].text(0.05, 0.95, f'r = {correlation:.3f}', transform=axes[1, 0].transAxes,
                       bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

        # 5. Degree distribution (log-log plot)
        degrees = [d for n, d in self.graph.degree()]
        degree_counts = Counter(degrees)
        x = list(degree_counts.keys())
        y = list(degree_counts.values())

        axes[1, 1].loglog(x, y, 'bo-', alpha=0.7, markersize=4)
        axes[1, 1].set_xlabel('Degree (log scale)')
        axes[1, 1].set_ylabel('Frequency (log scale)')
        axes[1, 1].set_title('Degree Distribution (Log-Log)')
        axes[1, 1].grid(True, alpha=0.3)

        # 6. Quality assessment summary
        if 'quality_distribution' in self.validation_results.get('precision_recall_metrics', {}):
            quality_data = self.validation_results['precision_recall_metrics']['quality_distribution']
            categories = ['High Quality', 'Medium Quality', 'Low Quality']
            values = [
                quality_data['high_quality'],
                quality_data['medium_quality'] - quality_data['high_quality'],
                quality_data['total_associations'] - quality_data['medium_quality']
            ]

            colors = ['darkgreen', 'orange', 'lightcoral']
            axes[1, 2].pie(values, labels=categories, colors=colors, autopct='%1.1f%%', startangle=90)
            axes[1, 2].set_title('Association Quality Distribution')

        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, 'validation_metrics_summary.png'), dpi=300, bbox_inches='tight')
        print(f"Validation visualizations saved to {output_dir}/validation_metrics_summary.png")

        return fig


def main():
    """Main function to run the advanced validation analysis"""

    # Set up paths
    base_dir = os.path.dirname(__file__)
    graph_path = os.path.join(base_dir, 'output', 'multi_gene_disease_graph.graphml')
    data_dir = os.path.join(base_dir, 'data')
    output_dir = os.path.join(base_dir, 'output')

    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)

    # Check if required files exist
    if not os.path.exists(graph_path):
        print(f"Error: Graph file not found at {graph_path}")
        print("Please run multi_gene_graph.py first to generate the graph.")
        return

    if not os.path.exists(data_dir):
        print(f"Error: Data directory not found at {data_dir}")
        return

    # Initialize validator
    print("Initializing Gene-Disease Precision Validator...")
    validator = GeneDiseasePrecisionValidator(graph_path, data_dir)

    # Generate comprehensive validation report
    results = validator.generate_comprehensive_report(output_dir)

    # Print summary
    print("\n" + "="*60)
    print("VALIDATION SUMMARY")
    print("="*60)

    if 'quality_metrics' in results:
        qm = results['quality_metrics']
        print(f"High-confidence associations: {qm.get('high_confidence_associations', 0)} ({qm.get('high_confidence_percentage', 0):.1f}%)")

    if 'topology_metrics' in results:
        tm = results['topology_metrics']
        print(f"Network density: {tm.get('density', 0):.4f}")
        print(f"Average clustering: {tm.get('avg_clustering', 0):.4f}")
        print(f"Scale-free network: {'Yes' if tm.get('is_scale_free', False) else 'No'}")

    if 'community_metrics' in results:
        cm = results['community_metrics']
        print(f"Modularity: {cm.get('modularity', 0):.4f}")
        print(f"Community purity: {cm.get('avg_community_purity', 0):.3f}")

    print("="*60)
    print("Validation analysis complete!")


if __name__ == "__main__":
    main()
