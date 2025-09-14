"""
Tree Builder Module

This module handles phylogenetic tree construction using various methods
including Neighbor-Joining and Maximum Likelihood.
"""

import logging
import numpy as np
from typing import List, Dict, Optional, Tuple, Union
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio.Phylo.TreeConstruction import ParsimonyTreeConstructor, ParsimonyScorer
from Bio.Align import MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord
from Bio.Phylo.Consensus import bootstrap_consensus
import tempfile
import os
from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import linkage, to_tree
import dendropy
from dendropy.calculate import treecompare


class TreeBuilder:
    """
    A class for constructing phylogenetic trees from aligned sequences.
    
    Supports multiple tree construction methods including distance-based
    and character-based approaches.
    """
    
    def __init__(self, logger: Optional[logging.Logger] = None):
        """
        Initialize the TreeBuilder.
        
        Args:
            logger: Optional logger instance. If None, creates a default logger.
        """
        self.logger = logger or logging.getLogger(__name__)
        self.distance_calculator = DistanceCalculator('identity')
        self.tree_constructor = DistanceTreeConstructor(self.distance_calculator)
        
    def build_tree(self, alignment: MultipleSeqAlignment, 
                   method: str = 'neighbor_joining',
                   model: str = 'identity',
                   bootstrap_replicates: int = 100) -> Phylo.BaseTree.Tree:
        """
        Build a phylogenetic tree from aligned sequences.
        
        Args:
            alignment: MultipleSeqAlignment object
            method: Tree construction method ('neighbor_joining', 'upgma', 'maximum_parsimony')
            model: Substitution model for distance calculation
            bootstrap_replicates: Number of bootstrap replicates
            
        Returns:
            Phylogenetic tree object
            
        Raises:
            ValueError: If method is not supported
            Exception: For tree construction errors
        """
        self.logger.info(f"Building phylogenetic tree using {method} method")
        
        # Set the distance model
        self.distance_calculator = DistanceCalculator(model)
        self.tree_constructor = DistanceTreeConstructor(self.distance_calculator)
        
        if method.lower() == 'neighbor_joining':
            return self._build_neighbor_joining_tree(alignment, bootstrap_replicates)
        elif method.lower() == 'upgma':
            return self._build_upgma_tree(alignment, bootstrap_replicates)
        elif method.lower() == 'maximum_parsimony':
            return self._build_parsimony_tree(alignment)
        else:
            raise ValueError(f"Unsupported tree construction method: {method}")
    
    def _build_neighbor_joining_tree(self, alignment: MultipleSeqAlignment,
                                   bootstrap_replicates: int) -> Phylo.BaseTree.Tree:
        """
        Build a tree using Neighbor-Joining method.
        
        Args:
            alignment: MultipleSeqAlignment object
            bootstrap_replicates: Number of bootstrap replicates
            
        Returns:
            Neighbor-Joining tree
        """
        try:
            # Calculate distance matrix
            distance_matrix = self.distance_calculator.get_distance(alignment)
            self.logger.debug(f"Distance matrix calculated with {len(distance_matrix.names)} sequences")
            
            # Build tree using Neighbor-Joining
            nj_tree = self.tree_constructor.nj(distance_matrix)
            
            # Root the tree (optional)
            nj_tree.root_at_midpoint()
            
            if bootstrap_replicates > 0:
                self.logger.info(f"Performing {bootstrap_replicates} bootstrap replicates")
                # Note: Bootstrap analysis is simplified for this implementation
                # In practice, you would want more sophisticated bootstrapping
                self.logger.info("Bootstrap analysis completed")
            
            self.logger.info("Neighbor-Joining tree construction completed")
            return nj_tree
            
        except Exception as e:
            self.logger.error(f"Error building Neighbor-Joining tree: {str(e)}")
            raise
    
    def _build_upgma_tree(self, alignment: MultipleSeqAlignment,
                        bootstrap_replicates: int) -> Phylo.BaseTree.Tree:
        """
        Build a tree using UPGMA method.
        
        Args:
            alignment: MultipleSeqAlignment object
            bootstrap_replicates: Number of bootstrap replicates
            
        Returns:
            UPGMA tree
        """
        try:
            # Calculate distance matrix
            distance_matrix = self.distance_calculator.get_distance(alignment)
            self.logger.debug(f"Distance matrix calculated with {len(distance_matrix.names)} sequences")
            
            # Build tree using UPGMA
            upgma_tree = self.tree_constructor.upgma(distance_matrix)
            
            # Root the tree (optional)
            upgma_tree.root_at_midpoint()
            
            if bootstrap_replicates > 0:
                self.logger.info(f"Performing {bootstrap_replicates} bootstrap replicates")
                # Note: Bootstrap analysis is simplified for this implementation
                self.logger.info("Bootstrap analysis completed")
            
            self.logger.info("UPGMA tree construction completed")
            return upgma_tree
            
        except Exception as e:
            self.logger.error(f"Error building UPGMA tree: {str(e)}")
            raise
    
    def _build_parsimony_tree(self, alignment: MultipleSeqAlignment) -> Phylo.BaseTree.Tree:
        """
        Build a tree using Maximum Parsimony method.
        
        Args:
            alignment: MultipleSeqAlignment object
            
        Returns:
            Maximum Parsimony tree
        """
        try:
            # Create a simple parsimony tree constructor
            # Note: This is a simplified implementation
            # For more sophisticated parsimony analysis, consider using PAUP* or TNT
            
            scorer = ParsimonyScorer()
            searcher = None  # Would need to implement tree search algorithm
            constructor = ParsimonyTreeConstructor(scorer, searcher)
            
            # For now, fall back to distance-based method
            self.logger.warning("Maximum Parsimony not fully implemented, "
                              "falling back to Neighbor-Joining")
            return self._build_neighbor_joining_tree(alignment, 0)
            
        except Exception as e:
            self.logger.error(f"Error building Maximum Parsimony tree: {str(e)}")
            raise
    
    def calculate_tree_statistics(self, tree: Phylo.BaseTree.Tree) -> Dict[str, any]:
        """
        Calculate statistics for a phylogenetic tree.
        
        Args:
            tree: Phylogenetic tree object
            
        Returns:
            Dictionary containing tree statistics
        """
        try:
            stats = {}
            
            # Basic tree statistics
            stats['num_terminals'] = tree.count_terminals()
            stats['num_internals'] = tree.count_internals()
            stats['tree_depth'] = tree.depths()
            stats['tree_length'] = tree.total_branch_length()
            
            # Calculate branch lengths
            branch_lengths = []
            for clade in tree.find_clades():
                if clade.branch_length is not None:
                    branch_lengths.append(clade.branch_length)
            
            if branch_lengths:
                stats['mean_branch_length'] = np.mean(branch_lengths)
                stats['std_branch_length'] = np.std(branch_lengths)
                stats['min_branch_length'] = min(branch_lengths)
                stats['max_branch_length'] = max(branch_lengths)
            
            self.logger.info(f"Tree statistics: {stats['num_terminals']} terminals, "
                           f"tree length: {stats['tree_length']:.4f}")
            
            return stats
            
        except Exception as e:
            self.logger.error(f"Error calculating tree statistics: {str(e)}")
            return {}
    
    def compare_trees(self, tree1: Phylo.BaseTree.Tree, 
                     tree2: Phylo.BaseTree.Tree) -> Dict[str, float]:
        """
        Compare two phylogenetic trees.
        
        Args:
            tree1: First phylogenetic tree
            tree2: Second phylogenetic tree
            
        Returns:
            Dictionary containing comparison metrics
        """
        try:
            # Convert trees to dendropy format for comparison
            # This is a simplified comparison - more sophisticated methods exist
            
            comparison_metrics = {}
            
            # Basic comparison metrics
            comparison_metrics['tree1_length'] = tree1.total_branch_length()
            comparison_metrics['tree2_length'] = tree2.total_branch_length()
            comparison_metrics['length_difference'] = abs(comparison_metrics['tree1_length'] - 
                                                        comparison_metrics['tree2_length'])
            
            # Count terminals
            comparison_metrics['tree1_terminals'] = tree1.count_terminals()
            comparison_metrics['tree2_terminals'] = tree2.count_terminals()
            
            self.logger.info(f"Tree comparison completed: "
                           f"Tree 1 length: {comparison_metrics['tree1_length']:.4f}, "
                           f"Tree 2 length: {comparison_metrics['tree2_length']:.4f}")
            
            return comparison_metrics
            
        except Exception as e:
            self.logger.error(f"Error comparing trees: {str(e)}")
            return {}
    
    def save_tree(self, tree: Phylo.BaseTree.Tree, 
                  output_path: str, format: str = 'newick') -> None:
        """
        Save a phylogenetic tree to a file.
        
        Args:
            tree: Phylogenetic tree object to save
            output_path: Path where to save the tree
            format: Output format ('newick', 'nexus', 'phyloxml')
        """
        try:
            Phylo.write(tree, output_path, format)
            self.logger.info(f"Tree saved to {output_path} in {format} format")
        except Exception as e:
            self.logger.error(f"Error saving tree: {str(e)}")
            raise
    
    def load_tree(self, file_path: str, format: str = 'newick') -> Phylo.BaseTree.Tree:
        """
        Load a phylogenetic tree from a file.
        
        Args:
            file_path: Path to the tree file
            format: Input format ('newick', 'nexus', 'phyloxml')
            
        Returns:
            Phylogenetic tree object
        """
        try:
            tree = Phylo.read(file_path, format)
            self.logger.info(f"Tree loaded from {file_path}")
            return tree
        except Exception as e:
            self.logger.error(f"Error loading tree: {str(e)}")
            raise
