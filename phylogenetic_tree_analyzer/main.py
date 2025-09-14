"""
Main Module for Phylogenetic Tree Analyzer

This module provides the command-line interface and orchestrates the
phylogenetic analysis workflow.
"""

import argparse
import logging
import sys
import os
from pathlib import Path
from typing import Optional, List
import time

from .sequence_parser import SequenceParser
from .alignment_manager import AlignmentManager
from .tree_builder import TreeBuilder
from .tree_visualizer import TreeVisualizer


def setup_logging(log_level: str = 'INFO', log_file: Optional[str] = None) -> logging.Logger:
    """
    Set up logging configuration.
    
    Args:
        log_level: Logging level (DEBUG, INFO, WARNING, ERROR, CRITICAL)
        log_file: Optional log file path
        
    Returns:
        Configured logger instance
    """
    # Create logger
    logger = logging.getLogger('phylo_analyzer')
    logger.setLevel(getattr(logging, log_level.upper()))
    
    # Clear any existing handlers
    logger.handlers.clear()
    
    # Create formatter
    formatter = logging.Formatter(
        '%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )
    
    # Console handler
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(getattr(logging, log_level.upper()))
    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)
    
    # File handler (if specified)
    if log_file:
        file_handler = logging.FileHandler(log_file)
        file_handler.setLevel(logging.DEBUG)  # Always log everything to file
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)
    
    return logger


def create_output_directories(output_dir: str, logger: logging.Logger) -> None:
    """
    Create output directories if they don't exist.
    
    Args:
        output_dir: Base output directory path
        logger: Logger instance
    """
    directories = [
        output_dir,
        os.path.join(output_dir, 'alignments'),
        os.path.join(output_dir, 'trees'),
        os.path.join(output_dir, 'plots'),
        os.path.join(output_dir, 'logs')
    ]
    
    for directory in directories:
        os.makedirs(directory, exist_ok=True)
        logger.debug(f"Created directory: {directory}")


def run_phylogenetic_analysis(input_file: str,
                             output_dir: str,
                             alignment_method: str = 'muscle',
                             tree_method: str = 'neighbor_joining',
                             substitution_model: str = 'identity',
                             bootstrap_replicates: int = 100,
                             min_sequence_length: int = 50,
                             log_level: str = 'INFO') -> None:
    """
    Run the complete phylogenetic analysis workflow.
    
    Args:
        input_file: Path to input sequence file
        output_dir: Output directory for results
        alignment_method: Alignment method ('muscle', 'clustalw', 'simple')
        tree_method: Tree construction method ('neighbor_joining', 'upgma', 'maximum_parsimony')
        substitution_model: Substitution model for distance calculation
        bootstrap_replicates: Number of bootstrap replicates
        min_sequence_length: Minimum sequence length filter
        log_level: Logging level
    """
    # Set up logging
    log_file = os.path.join(output_dir, 'logs', 'phylo_analysis.log')
    logger = setup_logging(log_level, log_file)
    
    try:
        logger.info("=" * 60)
        logger.info("PHYLOGENETIC TREE ANALYZER")
        logger.info("=" * 60)
        logger.info(f"Input file: {input_file}")
        logger.info(f"Output directory: {output_dir}")
        logger.info(f"Alignment method: {alignment_method}")
        logger.info(f"Tree method: {tree_method}")
        logger.info(f"Substitution model: {substitution_model}")
        logger.info(f"Bootstrap replicates: {bootstrap_replicates}")
        logger.info(f"Minimum sequence length: {min_sequence_length}")
        
        start_time = time.time()
        
        # Create output directories
        create_output_directories(output_dir, logger)
        
        # Step 1: Parse sequences
        logger.info("\n" + "-" * 40)
        logger.info("STEP 1: PARSING SEQUENCES")
        logger.info("-" * 40)
        
        parser = SequenceParser(logger)
        sequences = parser.parse_file(input_file)
        
        if not sequences:
            raise ValueError("No sequences found in input file")
        
        # Validate sequences
        sequences = parser.validate_sequences(sequences, min_sequence_length)
        
        if len(sequences) < 2:
            raise ValueError("At least 2 valid sequences are required for phylogenetic analysis")
        
        # Get sequence information
        seq_info = parser.get_sequence_info(sequences)
        logger.info(f"Sequences parsed successfully: {seq_info}")
        
        # Step 2: Multiple sequence alignment
        logger.info("\n" + "-" * 40)
        logger.info("STEP 2: MULTIPLE SEQUENCE ALIGNMENT")
        logger.info("-" * 40)
        
        aligner = AlignmentManager(logger)
        
        # Try the requested alignment method, fall back to simple if not available
        try:
            alignment = aligner.multiple_sequence_align(sequences, alignment_method)
            logger.info(f"Successfully used {alignment_method} for alignment")
        except Exception as e:
            logger.warning(f"{alignment_method} alignment failed: {str(e)}")
            logger.info("Falling back to simple alignment method")
            alignment = aligner.multiple_sequence_align(sequences, 'simple')
        
        # Save alignment
        alignment_output = os.path.join(output_dir, 'alignments', 'aligned_sequences.fasta')
        aligner.save_alignment(alignment, alignment_output)
        
        # Calculate alignment statistics
        alignment_stats = aligner.calculate_alignment_statistics(alignment)
        logger.info(f"Alignment completed: {alignment_stats}")
        
        # Step 3: Phylogenetic tree construction
        logger.info("\n" + "-" * 40)
        logger.info("STEP 3: PHYLOGENETIC TREE CONSTRUCTION")
        logger.info("-" * 40)
        
        tree_builder = TreeBuilder(logger)
        tree = tree_builder.build_tree(alignment, tree_method, substitution_model, bootstrap_replicates)
        
        # Save tree
        tree_output = os.path.join(output_dir, 'trees', f'phylogenetic_tree_{tree_method}.newick')
        tree_builder.save_tree(tree, tree_output)
        
        # Calculate tree statistics
        tree_stats = tree_builder.calculate_tree_statistics(tree)
        logger.info(f"Tree construction completed: {tree_stats}")
        
        # Step 4: Visualization
        logger.info("\n" + "-" * 40)
        logger.info("STEP 4: TREE VISUALIZATION")
        logger.info("-" * 40)
        
        visualizer = TreeVisualizer(logger)
        
        # Create tree plot
        plot_output = os.path.join(output_dir, 'plots', f'phylogenetic_tree_{tree_method}.png')
        fig = visualizer.plot_tree(
            tree,
            title=f"Phylogenetic Tree ({tree_method.replace('_', ' ').title()})",
            output_path=plot_output,
            figsize=(14, 10),
            show_branch_lengths=True,
            show_bootstrap_values=True
        )
        
        # Create summary plot
        summary_output = os.path.join(output_dir, 'plots', 'analysis_summary.png')
        visualizer.create_summary_plot(alignment_stats, tree_stats, summary_output)
        
        # Step 5: Generate report
        logger.info("\n" + "-" * 40)
        logger.info("STEP 5: GENERATING REPORT")
        logger.info("-" * 40)
        
        report_path = os.path.join(output_dir, 'analysis_report.txt')
        generate_analysis_report(
            input_file, output_dir, alignment_stats, tree_stats, 
            alignment_method, tree_method, report_path, logger
        )
        
        # Calculate total runtime
        end_time = time.time()
        runtime = end_time - start_time
        
        logger.info("\n" + "=" * 60)
        logger.info("ANALYSIS COMPLETED SUCCESSFULLY")
        logger.info("=" * 60)
        logger.info(f"Total runtime: {runtime:.2f} seconds")
        logger.info(f"Results saved to: {output_dir}")
        logger.info("Files generated:")
        logger.info(f"  - Alignment: {alignment_output}")
        logger.info(f"  - Tree: {tree_output}")
        logger.info(f"  - Tree plot: {plot_output}")
        logger.info(f"  - Summary plot: {summary_output}")
        logger.info(f"  - Report: {report_path}")
        logger.info(f"  - Log file: {log_file}")
        
    except Exception as e:
        logger.error(f"Analysis failed: {str(e)}")
        logger.exception("Full traceback:")
        raise


def generate_analysis_report(input_file: str, output_dir: str,
                           alignment_stats: dict, tree_stats: dict,
                           alignment_method: str, tree_method: str,
                           report_path: str, logger: logging.Logger) -> None:
    """
    Generate a text report summarizing the analysis results.
    
    Args:
        input_file: Input sequence file path
        output_dir: Output directory
        alignment_stats: Alignment statistics
        tree_stats: Tree statistics
        alignment_method: Alignment method used
        tree_method: Tree construction method used
        report_path: Path to save the report
        logger: Logger instance
    """
    try:
        with open(report_path, 'w') as f:
            f.write("PHYLOGENETIC TREE ANALYSIS REPORT\n")
            f.write("=" * 50 + "\n\n")
            
            f.write(f"Analysis Date: {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"Input File: {input_file}\n")
            f.write(f"Output Directory: {output_dir}\n\n")
            
            f.write("ANALYSIS PARAMETERS\n")
            f.write("-" * 20 + "\n")
            f.write(f"Alignment Method: {alignment_method}\n")
            f.write(f"Tree Construction Method: {tree_method}\n\n")
            
            f.write("SEQUENCE ALIGNMENT STATISTICS\n")
            f.write("-" * 30 + "\n")
            for key, value in alignment_stats.items():
                if isinstance(value, float):
                    f.write(f"{key.replace('_', ' ').title()}: {value:.3f}\n")
                else:
                    f.write(f"{key.replace('_', ' ').title()}: {value}\n")
            f.write("\n")
            
            f.write("PHYLOGENETIC TREE STATISTICS\n")
            f.write("-" * 30 + "\n")
            for key, value in tree_stats.items():
                if isinstance(value, float):
                    f.write(f"{key.replace('_', ' ').title()}: {value:.4f}\n")
                else:
                    f.write(f"{key.replace('_', ' ').title()}: {value}\n")
            f.write("\n")
            
            f.write("OUTPUT FILES\n")
            f.write("-" * 12 + "\n")
            f.write(f"Aligned sequences: {os.path.join(output_dir, 'alignments', 'aligned_sequences.fasta')}\n")
            f.write(f"Phylogenetic tree: {os.path.join(output_dir, 'trees', f'phylogenetic_tree_{tree_method}.newick')}\n")
            f.write(f"Tree visualization: {os.path.join(output_dir, 'plots', f'phylogenetic_tree_{tree_method}.png')}\n")
            f.write(f"Summary plot: {os.path.join(output_dir, 'plots', 'analysis_summary.png')}\n")
            f.write(f"Log file: {os.path.join(output_dir, 'logs', 'phylo_analysis.log')}\n")
        
        logger.info(f"Analysis report saved to {report_path}")
        
    except Exception as e:
        logger.error(f"Error generating report: {str(e)}")
        raise


def main():
    """Main entry point for the command-line interface."""
    parser = argparse.ArgumentParser(
        description="Phylogenetic Tree Analyzer - Analyze evolutionary relationships among species",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic analysis with default parameters
  python -m phylogenetic_tree_analyzer data/sample_sequences.fasta -o output/

  # Custom analysis with specific methods
  python -m phylogenetic_tree_analyzer data/sample_sequences.fasta -o output/ \\
    --alignment-method muscle --tree-method neighbor_joining \\
    --bootstrap-replicates 1000 --log-level DEBUG

  # Analyze GenBank file
  python -m phylogenetic_tree_analyzer data/sample_sequences.gb -o output/ \\
    --alignment-method clustalw --tree-method upgma
        """
    )
    
    # Required arguments
    parser.add_argument('input_file',
                       help='Path to input sequence file (FASTA or GenBank format)')
    parser.add_argument('-o', '--output-dir',
                       required=True,
                       help='Output directory for results')
    
    # Analysis parameters
    parser.add_argument('--alignment-method',
                       choices=['muscle', 'clustalw', 'simple'],
                       default='muscle',
                       help='Multiple sequence alignment method (default: muscle)')
    
    parser.add_argument('--tree-method',
                       choices=['neighbor_joining', 'upgma', 'maximum_parsimony'],
                       default='neighbor_joining',
                       help='Phylogenetic tree construction method (default: neighbor_joining)')
    
    parser.add_argument('--substitution-model',
                       choices=['identity', 'blastn', 'trans'],
                       default='identity',
                       help='Substitution model for distance calculation (default: identity)')
    
    parser.add_argument('--bootstrap-replicates',
                       type=int,
                       default=100,
                       help='Number of bootstrap replicates (default: 100)')
    
    parser.add_argument('--min-sequence-length',
                       type=int,
                       default=50,
                       help='Minimum sequence length for analysis (default: 50)')
    
    # Output options
    parser.add_argument('--log-level',
                       choices=['DEBUG', 'INFO', 'WARNING', 'ERROR'],
                       default='INFO',
                       help='Logging level (default: INFO)')
    
    parser.add_argument('--version',
                       action='version',
                       version='%(prog)s 1.0.0')
    
    # Parse arguments
    args = parser.parse_args()
    
    # Validate input file
    if not os.path.exists(args.input_file):
        print(f"Error: Input file '{args.input_file}' not found.")
        sys.exit(1)
    
    # Run analysis
    try:
        run_phylogenetic_analysis(
            input_file=args.input_file,
            output_dir=args.output_dir,
            alignment_method=args.alignment_method,
            tree_method=args.tree_method,
            substitution_model=args.substitution_model,
            bootstrap_replicates=args.bootstrap_replicates,
            min_sequence_length=args.min_sequence_length,
            log_level=args.log_level
        )
        
        print(f"\nAnalysis completed successfully! Results saved to: {args.output_dir}")
        
    except KeyboardInterrupt:
        print("\nAnalysis interrupted by user.")
        sys.exit(1)
    except Exception as e:
        print(f"\nAnalysis failed: {str(e)}")
        sys.exit(1)


if __name__ == '__main__':
    main()
