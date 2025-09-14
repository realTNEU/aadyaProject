#!/usr/bin/env python3
"""
Example usage script for Phylogenetic Tree Analyzer.

This script demonstrates how to use the package programmatically
without the command-line interface.
"""

import os
import sys
from pathlib import Path

# Add the current directory to Python path to import the package
sys.path.insert(0, str(Path(__file__).parent))

from phylogenetic_tree_analyzer import (
    SequenceParser, AlignmentManager, TreeBuilder, TreeVisualizer
)


def run_example_analysis(input_file: str, output_dir: str = "example_output"):
    """
    Run a complete phylogenetic analysis example.
    
    Args:
        input_file: Path to input FASTA file
        output_dir: Output directory for results
    """
    print("=" * 60)
    print("PHYLOGENETIC TREE ANALYZER - EXAMPLE USAGE")
    print("=" * 60)
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    try:
        # Step 1: Parse sequences
        print("\nStep 1: Parsing sequences...")
        parser = SequenceParser()
        sequences = parser.parse_file(input_file)
        
        if not sequences:
            print("No sequences found in input file!")
            return
        
        print(f"Parsed {len(sequences)} sequences")
        
        # Validate sequences
        sequences = parser.validate_sequences(sequences, min_length=50)
        print(f"Validated {len(sequences)} sequences")
        
        # Get sequence information
        seq_info = parser.get_sequence_info(sequences)
        print(f"Sequence info: {seq_info}")
        
        # Step 2: Multiple sequence alignment
        print("\nStep 2: Multiple sequence alignment...")
        aligner = AlignmentManager()
        
        # Try MUSCLE first, fall back to simple alignment if not available
        try:
            alignment = aligner.multiple_sequence_align(sequences, method='muscle')
            print("Used MUSCLE for alignment")
        except:
            print("MUSCLE not available, using simple alignment")
            alignment = aligner.multiple_sequence_align(sequences, method='simple')
        
        # Save alignment
        alignment_output = os.path.join(output_dir, "aligned_sequences.fasta")
        aligner.save_alignment(alignment, alignment_output)
        print(f"Alignment saved to: {alignment_output}")
        
        # Calculate alignment statistics
        alignment_stats = aligner.calculate_alignment_statistics(alignment)
        print(f"Alignment statistics: {alignment_stats}")
        
        # Step 3: Phylogenetic tree construction
        print("\nStep 3: Phylogenetic tree construction...")
        tree_builder = TreeBuilder()
        tree = tree_builder.build_tree(alignment, method='neighbor_joining')
        
        # Save tree
        tree_output = os.path.join(output_dir, "phylogenetic_tree.newick")
        tree_builder.save_tree(tree, tree_output)
        print(f"Tree saved to: {tree_output}")
        
        # Calculate tree statistics
        tree_stats = tree_builder.calculate_tree_statistics(tree)
        print(f"Tree statistics: {tree_stats}")
        
        # Step 4: Visualization
        print("\nStep 4: Tree visualization...")
        visualizer = TreeVisualizer()
        
        # Create tree plot
        plot_output = os.path.join(output_dir, "phylogenetic_tree.png")
        fig = visualizer.plot_tree(
            tree,
            title="Example Phylogenetic Tree",
            output_path=plot_output,
            figsize=(12, 8)
        )
        print(f"Tree plot saved to: {plot_output}")
        
        # Create summary plot
        summary_output = os.path.join(output_dir, "analysis_summary.png")
        visualizer.create_summary_plot(alignment_stats, tree_stats, summary_output)
        print(f"Summary plot saved to: {summary_output}")
        
        print("\n" + "=" * 60)
        print("ANALYSIS COMPLETED SUCCESSFULLY!")
        print("=" * 60)
        print(f"Results saved to: {output_dir}")
        print("\nGenerated files:")
        print(f"  - Alignment: {alignment_output}")
        print(f"  - Tree: {tree_output}")
        print(f"  - Tree plot: {plot_output}")
        print(f"  - Summary plot: {summary_output}")
        
    except Exception as e:
        print(f"\nError during analysis: {str(e)}")
        import traceback
        traceback.print_exc()


if __name__ == "__main__":
    # Check if input file is provided
    if len(sys.argv) < 2:
        print("Usage: python example_usage.py <input_fasta_file> [output_directory]")
        print("\nExample:")
        print("  python example_usage.py data/sample_sequences.fasta")
        print("  python example_usage.py data/sample_sequences.fasta my_output")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_dir = sys.argv[2] if len(sys.argv) > 2 else "example_output"
    
    # Check if input file exists
    if not os.path.exists(input_file):
        print(f"Error: Input file '{input_file}' not found.")
        sys.exit(1)
    
    # Run the analysis
    run_example_analysis(input_file, output_dir)
